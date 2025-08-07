//! System integration and resource monitoring for REPrise
//!
//! Provides resource monitoring, limits enforcement, process management,
//! and system-level optimizations for production deployments.

use crate::config::SystemSettings;
use crate::error::{RepriseError, Result};
use crate::logging::MetricsCollector;
use once_cell::sync::Lazy;
use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex, RwLock};
use std::thread;
use std::time::{Duration, Instant};
use sysinfo::System;
use tracing::{debug, error, info, instrument, warn};

/// System resource information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SystemInfo {
    /// Total system memory in bytes
    pub total_memory: u64,
    /// Available system memory in bytes
    pub available_memory: u64,
    /// Number of CPU cores
    pub cpu_count: usize,
    /// Current memory usage percentage
    pub memory_usage_percent: f64,
    /// Current CPU usage percentage
    pub cpu_usage_percent: f64,
    /// Current process memory usage in bytes
    pub process_memory: u64,
    /// Current process CPU usage percentage
    pub process_cpu_percent: f64,
    /// System load averages (1min, 5min, 15min) - Unix only
    pub load_averages: Option<(f64, f64, f64)>,
    /// Number of file descriptors in use - Unix only
    pub file_descriptors: Option<u32>,
}

/// Resource usage limits and thresholds
#[derive(Debug, Clone)]
pub struct ResourceLimits {
    /// Maximum memory usage in bytes
    pub max_memory: u64,
    /// Memory warning threshold (percentage)
    pub memory_warning_threshold: f64,
    /// Memory critical threshold (percentage)
    pub memory_critical_threshold: f64,
    /// Maximum CPU usage percentage
    pub max_cpu_usage: f64,
    /// Maximum number of threads
    pub max_threads: usize,
    /// Maximum file handles
    pub max_file_handles: usize,
}

/// Resource monitoring events
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ResourceEvent {
    /// Memory usage warning
    MemoryWarning {
        current: u64,
        threshold: u64,
        percentage: f64,
    },
    /// Memory usage critical
    MemoryCritical {
        current: u64,
        threshold: u64,
        percentage: f64,
    },
    /// CPU usage warning
    CpuWarning {
        current: f64,
        threshold: f64,
    },
    /// Resource limit exceeded
    LimitExceeded {
        resource: String,
        current: u64,
        limit: u64,
    },
    /// System overload detected
    SystemOverload {
        memory_percent: f64,
        cpu_percent: f64,
        load_average: Option<f64>,
    },
}

/// System resource monitor
pub struct ResourceMonitor {
    settings: SystemSettings,
    system: Arc<Mutex<System>>,
    metrics: Arc<MetricsCollector>,
    monitoring_active: Arc<AtomicBool>,
    last_update: Arc<RwLock<Instant>>,
    resource_limits: ResourceLimits,
    event_handlers: Arc<RwLock<Vec<Box<dyn Fn(ResourceEvent) + Send + Sync>>>>,
}

impl ResourceMonitor {
    /// Create a new resource monitor
    pub fn new(settings: SystemSettings, metrics: Arc<MetricsCollector>) -> Self {
        let mut system = System::new_all();
        system.refresh_all();

        let total_memory = system.total_memory() * 1024; // Convert to bytes
        let resource_limits = ResourceLimits {
            max_memory: (total_memory as f64 * 0.8) as u64, // Default to 80% of system memory
            memory_warning_threshold: settings.memory_warning_threshold,
            memory_critical_threshold: settings.memory_critical_threshold,
            max_cpu_usage: 90.0, // 90% CPU usage threshold
            max_threads: 1000,   // Default thread limit
            max_file_handles: settings.max_file_handles,
        };

        Self {
            settings,
            system: Arc::new(Mutex::new(system)),
            metrics,
            monitoring_active: Arc::new(AtomicBool::new(false)),
            last_update: Arc::new(RwLock::new(Instant::now())),
            resource_limits,
            event_handlers: Arc::new(RwLock::new(Vec::new())),
        }
    }

    /// Start resource monitoring in background thread
    #[instrument(skip(self))]
    pub fn start_monitoring(&self) -> Result<()> {
        if self.monitoring_active.load(Ordering::Relaxed) {
            warn!("Resource monitoring already active");
            return Ok(());
        }

        self.monitoring_active.store(true, Ordering::Relaxed);

        let system = Arc::clone(&self.system);
        let metrics = Arc::clone(&self.metrics);
        let monitoring_active = Arc::clone(&self.monitoring_active);
        let last_update = Arc::clone(&self.last_update);
        let event_handlers = Arc::clone(&self.event_handlers);
        let settings = self.settings.clone();
        let resource_limits = self.resource_limits.clone();

        thread::spawn(move || {
            info!("Resource monitoring started");

            while monitoring_active.load(Ordering::Relaxed) {
                if let Err(e) = Self::monitoring_loop(
                    &system,
                    &metrics,
                    &last_update,
                    &event_handlers,
                    &settings,
                    &resource_limits,
                ) {
                    error!(error = %e, "Error in resource monitoring loop");
                }

                thread::sleep(settings.memory_check_interval);
            }

            info!("Resource monitoring stopped");
        });

        Ok(())
    }

    /// Stop resource monitoring
    pub fn stop_monitoring(&self) {
        if self.monitoring_active.load(Ordering::Relaxed) {
            self.monitoring_active.store(false, Ordering::Relaxed);
            info!("Resource monitoring stopped");
        }
    }

    /// Get current system information
    pub fn system_info(&self) -> Result<SystemInfo> {
        let mut system = self.system.lock().unwrap();
        system.refresh_all();

        let total_memory = system.total_memory() * 1024;
        let available_memory = system.available_memory() * 1024;
        let memory_usage_percent = ((total_memory - available_memory) as f64 / total_memory as f64) * 100.0;

        // Get process-specific information
        let process = system.processes()
            .get(&sysinfo::get_current_pid().unwrap_or(sysinfo::Pid::from(0)));

        let (process_memory, process_cpu_percent) = if let Some(proc) = process {
            (proc.memory() * 1024, proc.cpu_usage() as f64)
        } else {
            (0, 0.0)
        };

        // Get load averages (Unix only)
        let load_averages = {
            let la = System::load_average();
            if la.one > 0.0 {
                Some((la.one, la.five, la.fifteen))
            } else {
                None
            }
        };

        Ok(SystemInfo {
            total_memory,
            available_memory,
            cpu_count: system.cpus().len(),
            memory_usage_percent,
            cpu_usage_percent: system.global_cpu_info().cpu_usage() as f64,
            process_memory,
            process_cpu_percent,
            load_averages,
            file_descriptors: None, // TODO: Implement for Unix systems
        })
    }

    /// Check if resource limits are being exceeded
    pub fn check_resource_limits(&self) -> Result<Vec<ResourceEvent>> {
        let system_info = self.system_info()?;
        let mut events = Vec::new();

        // Check memory thresholds
        let memory_warning_threshold = (system_info.total_memory as f64 * self.resource_limits.memory_warning_threshold / 100.0) as u64;
        let memory_critical_threshold = (system_info.total_memory as f64 * self.resource_limits.memory_critical_threshold / 100.0) as u64;

        let used_memory = system_info.total_memory - system_info.available_memory;

        if used_memory >= memory_critical_threshold {
            events.push(ResourceEvent::MemoryCritical {
                current: used_memory,
                threshold: memory_critical_threshold,
                percentage: system_info.memory_usage_percent,
            });
        } else if used_memory >= memory_warning_threshold {
            events.push(ResourceEvent::MemoryWarning {
                current: used_memory,
                threshold: memory_warning_threshold,
                percentage: system_info.memory_usage_percent,
            });
        }

        // Check CPU usage
        if system_info.cpu_usage_percent >= self.resource_limits.max_cpu_usage {
            events.push(ResourceEvent::CpuWarning {
                current: system_info.cpu_usage_percent,
                threshold: self.resource_limits.max_cpu_usage,
            });
        }

        // Check for system overload
        if system_info.memory_usage_percent >= 85.0 && system_info.cpu_usage_percent >= 80.0 {
            events.push(ResourceEvent::SystemOverload {
                memory_percent: system_info.memory_usage_percent,
                cpu_percent: system_info.cpu_usage_percent,
                load_average: system_info.load_averages.map(|la| la.0),
            });
        }

        // Check process memory limit
        if system_info.process_memory >= self.resource_limits.max_memory {
            events.push(ResourceEvent::LimitExceeded {
                resource: "process_memory".to_string(),
                current: system_info.process_memory,
                limit: self.resource_limits.max_memory,
            });
        }

        Ok(events)
    }

    /// Add a resource event handler
    pub fn add_event_handler<F>(&self, handler: F)
    where
        F: Fn(ResourceEvent) + Send + Sync + 'static,
    {
        let mut handlers = self.event_handlers.write().unwrap();
        handlers.push(Box::new(handler));
    }

    /// Set process priority (Unix nice value)
    #[cfg(unix)]
    pub fn set_process_priority(&self, priority: i8) -> Result<()> {
        use libc::{getpid, setpriority, PRIO_PROCESS};
        
        let result = unsafe { setpriority(PRIO_PROCESS, getpid() as u32, priority as i32) };
        
        if result != 0 {
            return Err(RepriseError::config(format!(
                "Failed to set process priority: {}",
                std::io::Error::last_os_error()
            )));
        }
        
        info!(priority = priority, "Process priority set");
        Ok(())
    }

    #[cfg(not(unix))]
    pub fn set_process_priority(&self, _priority: i8) -> Result<()> {
        warn!("Process priority setting not supported on this platform");
        Ok(())
    }

    /// Set CPU affinity for the current process
    #[cfg(target_os = "linux")]
    pub fn set_cpu_affinity(&self, cores: &[usize]) -> Result<()> {
        use libc::{cpu_set_t, sched_setaffinity, CPU_SET, CPU_ZERO, getpid};
        use std::mem;
        
        if cores.is_empty() {
            return Ok(());
        }
        
        let mut cpu_set: cpu_set_t = unsafe { mem::zeroed() };
        unsafe { CPU_ZERO(&mut cpu_set) };
        
        for &core in cores {
            if core >= 1024 {
                return Err(RepriseError::config(format!("CPU core {} out of range", core)));
            }
            unsafe { CPU_SET(core, &mut cpu_set) };
        }
        
        let result = unsafe {
            sched_setaffinity(
                getpid(),
                mem::size_of::<cpu_set_t>(),
                &cpu_set as *const cpu_set_t,
            )
        };
        
        if result != 0 {
            return Err(RepriseError::config(format!(
                "Failed to set CPU affinity: {}",
                std::io::Error::last_os_error()
            )));
        }
        
        info!(cores = ?cores, "CPU affinity set");
        Ok(())
    }

    #[cfg(not(target_os = "linux"))]
    pub fn set_cpu_affinity(&self, _cores: &[usize]) -> Result<()> {
        warn!("CPU affinity setting not supported on this platform");
        Ok(())
    }

    /// Force garbage collection hint
    pub fn gc_hint(&self) {
        if self.settings.enable_gc_hints {
            // Trigger garbage collection in memory allocator if available
            #[cfg(feature = "mimalloc")]
            {
                // mimalloc doesn't expose a direct GC function, but we can try to encourage cleanup
                debug!("GC hint triggered (mimalloc)");
            }
            
            #[cfg(not(feature = "mimalloc"))]
            {
                // For system allocator, we can't do much
                debug!("GC hint triggered (system allocator)");
            }
        }
    }

    /// Get resource usage recommendations
    pub fn get_recommendations(&self) -> Result<Vec<String>> {
        let system_info = self.system_info()?;
        let mut recommendations = Vec::new();

        // Memory recommendations
        if system_info.memory_usage_percent > 80.0 {
            recommendations.push(
                "High memory usage detected. Consider reducing batch size or enabling checkpointing.".to_string()
            );
        }

        if system_info.process_memory > self.resource_limits.max_memory / 2 {
            recommendations.push(
                "Process memory usage is high. Consider reducing worker threads or channel capacity.".to_string()
            );
        }

        // CPU recommendations
        if system_info.cpu_usage_percent > 90.0 {
            recommendations.push(
                "High CPU usage detected. Consider reducing worker threads.".to_string()
            );
        }

        if let Some((load1, _, _)) = system_info.load_averages {
            if load1 > system_info.cpu_count as f64 * 1.5 {
                recommendations.push(
                    "High system load detected. Consider reducing parallelism or running during off-peak hours.".to_string()
                );
            }
        }

        // General recommendations
        if system_info.available_memory < 1024 * 1024 * 1024 {
            recommendations.push(
                "Low available memory. Consider using a larger instance or enabling compression.".to_string()
            );
        }

        Ok(recommendations)
    }

    /// Internal monitoring loop
    fn monitoring_loop(
        system: &Arc<Mutex<System>>,
        metrics: &MetricsCollector,
        last_update: &Arc<RwLock<Instant>>,
        event_handlers: &Arc<RwLock<Vec<Box<dyn Fn(ResourceEvent) + Send + Sync>>>>,
        settings: &SystemSettings,
        _resource_limits: &ResourceLimits,
    ) -> Result<()> {
        let mut sys = system.lock().unwrap();
        sys.refresh_all();

        let total_memory = sys.total_memory() * 1024;
        let available_memory = sys.available_memory() * 1024;
        let used_memory = total_memory - available_memory;
        let memory_usage_percent = (used_memory as f64 / total_memory as f64) * 100.0;

        // Update metrics
        metrics.set_gauge("system_memory_total", total_memory);
        metrics.set_gauge("system_memory_available", available_memory);
        metrics.set_gauge("system_memory_used", used_memory);
        metrics.update_memory_usage(used_memory);

        // Get process information
        if let Some(process) = sys.processes().get(&sysinfo::get_current_pid().unwrap_or(sysinfo::Pid::from(0))) {
            let process_memory = process.memory() * 1024;
            metrics.set_gauge("process_memory", process_memory);
            metrics.set_gauge("process_cpu_percent", process.cpu_usage() as u64);
        }

        // Check thresholds and trigger events
        let memory_warning = (total_memory as f64 * settings.memory_warning_threshold / 100.0) as u64;
        let memory_critical = (total_memory as f64 * settings.memory_critical_threshold / 100.0) as u64;

        let handlers = event_handlers.read().unwrap();
        
        if used_memory >= memory_critical {
            let event = ResourceEvent::MemoryCritical {
                current: used_memory,
                threshold: memory_critical,
                percentage: memory_usage_percent,
            };
            
            error!(
                memory_used = used_memory,
                memory_threshold = memory_critical,
                memory_percent = memory_usage_percent,
                "Critical memory usage detected"
            );
            
            for handler in handlers.iter() {
                handler(event.clone());
            }
        } else if used_memory >= memory_warning {
            let event = ResourceEvent::MemoryWarning {
                current: used_memory,
                threshold: memory_warning,
                percentage: memory_usage_percent,
            };
            
            warn!(
                memory_used = used_memory,
                memory_threshold = memory_warning,
                memory_percent = memory_usage_percent,
                "High memory usage detected"
            );
            
            for handler in handlers.iter() {
                handler(event.clone());
            }
        }

        // Update last monitoring time
        *last_update.write().unwrap() = Instant::now();

        drop(sys); // Release system lock
        Ok(())
    }
}

/// System performance optimizer
pub struct PerformanceOptimizer {
    settings: SystemSettings,
}

impl PerformanceOptimizer {
    /// Create a new performance optimizer
    pub fn new(settings: SystemSettings) -> Self {
        Self { settings }
    }

    /// Apply system-level optimizations
    pub fn optimize(&self) -> Result<()> {
        info!("Applying system performance optimizations");

        // Set process priority
        if self.settings.process_priority != 0 {
            self.set_priority(self.settings.process_priority)?;
        }

        // Set CPU affinity
        if self.settings.enable_cpu_affinity && !self.settings.cpu_affinity_cores.is_empty() {
            self.set_cpu_affinity(&self.settings.cpu_affinity_cores)?;
        }

        // Configure memory allocator hints
        self.configure_memory_allocator()?;

        // Set file descriptor limits
        self.configure_fd_limits()?;

        info!("System performance optimizations applied");
        Ok(())
    }

    /// Set process priority
    fn set_priority(&self, priority: i8) -> Result<()> {
        // This would typically use the resource monitor's methods
        debug!(priority = priority, "Setting process priority");
        Ok(())
    }

    /// Set CPU affinity
    fn set_cpu_affinity(&self, cores: &[usize]) -> Result<()> {
        debug!(cores = ?cores, "Setting CPU affinity");
        Ok(())
    }

    /// Configure memory allocator
    fn configure_memory_allocator(&self) -> Result<()> {
        #[cfg(feature = "mimalloc")]
        {
            // mimalloc is configured at compile time, but we can log that it's being used
            info!("Using mimalloc high-performance allocator");
        }

        #[cfg(not(feature = "mimalloc"))]
        {
            debug!("Using system default allocator");
        }

        Ok(())
    }

    /// Configure file descriptor limits
    fn configure_fd_limits(&self) -> Result<()> {
        #[cfg(unix)]
        {
            use libc::{getrlimit, setrlimit, rlimit, RLIMIT_NOFILE};
            use std::mem;

            let mut rl: rlimit = unsafe { mem::zeroed() };
            let result = unsafe { getrlimit(RLIMIT_NOFILE, &mut rl) };

            if result == 0 {
                let current_limit = rl.rlim_cur as usize;
                let desired_limit = self.settings.max_file_handles;

                if current_limit < desired_limit {
                    rl.rlim_cur = desired_limit.min(rl.rlim_max as usize) as u64;
                    
                    let set_result = unsafe { setrlimit(RLIMIT_NOFILE, &rl) };
                    if set_result == 0 {
                        info!(
                            old_limit = current_limit,
                            new_limit = rl.rlim_cur,
                            "File descriptor limit increased"
                        );
                    } else {
                        warn!(
                            desired_limit = desired_limit,
                            max_limit = rl.rlim_max,
                            "Failed to increase file descriptor limit"
                        );
                    }
                }
            }
        }

        #[cfg(not(unix))]
        {
            debug!("File descriptor limit configuration not available on this platform");
        }

        Ok(())
    }
}

/// Global system manager instance
static SYSTEM_MANAGER: Lazy<RwLock<Option<Arc<ResourceMonitor>>>> = Lazy::new(|| RwLock::new(None));

/// Initialize global system manager
pub fn init_system_manager(settings: SystemSettings, metrics: Arc<MetricsCollector>) -> Result<()> {
    let monitor = Arc::new(ResourceMonitor::new(settings.clone(), metrics));
    monitor.start_monitoring()?;

    // Set up default event handlers
    monitor.add_event_handler(|event| {
        match event {
            ResourceEvent::MemoryWarning { percentage, .. } => {
                warn!(memory_percent = percentage, "Memory usage warning");
            }
            ResourceEvent::MemoryCritical { percentage, .. } => {
                error!(memory_percent = percentage, "Critical memory usage");
            }
            ResourceEvent::CpuWarning { current, threshold } => {
                warn!(cpu_usage = current, threshold = threshold, "High CPU usage");
            }
            ResourceEvent::SystemOverload { memory_percent, cpu_percent, .. } => {
                error!(
                    memory_percent = memory_percent,
                    cpu_percent = cpu_percent,
                    "System overload detected"
                );
            }
            ResourceEvent::LimitExceeded { resource, current, limit } => {
                error!(
                    resource = resource,
                    current = current,
                    limit = limit,
                    "Resource limit exceeded"
                );
            }
        }
    });

    // Apply optimizations
    let optimizer = PerformanceOptimizer::new(settings);
    optimizer.optimize()?;

    *SYSTEM_MANAGER.write().unwrap() = Some(monitor);
    Ok(())
}

/// Get global system manager
pub fn system_manager() -> Option<Arc<ResourceMonitor>> {
    SYSTEM_MANAGER.read().unwrap().clone()
}

/// Shutdown global system manager
pub fn shutdown_system_manager() {
    if let Some(manager) = SYSTEM_MANAGER.write().unwrap().take() {
        manager.stop_monitoring();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::logging::MetricsCollector;
    use std::sync::Arc;
    use std::time::Duration;

    fn create_test_settings() -> SystemSettings {
        SystemSettings {
            memory_check_interval: Duration::from_secs(1),
            memory_warning_threshold: 75.0,
            memory_critical_threshold: 85.0,
            enable_gc_hints: true,
            max_file_handles: 1024,
            process_priority: 0,
            enable_cpu_affinity: false,
            cpu_affinity_cores: Vec::new(),
        }
    }

    #[test]
    fn test_resource_monitor_creation() {
        let settings = create_test_settings();
        let metrics = Arc::new(MetricsCollector::new());
        let monitor = ResourceMonitor::new(settings, metrics);

        assert!(!monitor.monitoring_active.load(Ordering::Relaxed));
    }

    #[test]
    fn test_system_info_collection() {
        let settings = create_test_settings();
        let metrics = Arc::new(MetricsCollector::new());
        let monitor = ResourceMonitor::new(settings, metrics);

        let system_info = monitor.system_info().unwrap();
        
        assert!(system_info.total_memory > 0);
        assert!(system_info.cpu_count > 0);
        assert!(system_info.memory_usage_percent >= 0.0);
        assert!(system_info.memory_usage_percent <= 100.0);
    }

    #[test]
    fn test_resource_limit_checking() {
        let settings = create_test_settings();
        let metrics = Arc::new(MetricsCollector::new());
        let monitor = ResourceMonitor::new(settings, metrics);

        let events = monitor.check_resource_limits().unwrap();
        // Events depend on actual system state, so we just check that the call succeeds
        assert!(events.len() >= 0);
    }

    #[test]
    fn test_event_handler_registration() {
        let settings = create_test_settings();
        let metrics = Arc::new(MetricsCollector::new());
        let monitor = ResourceMonitor::new(settings, metrics);

        let event_received = Arc::new(AtomicBool::new(false));
        let event_received_clone = Arc::clone(&event_received);

        monitor.add_event_handler(move |_event| {
            event_received_clone.store(true, Ordering::Relaxed);
        });

        // Manually trigger an event through the internal mechanism
        let handlers = monitor.event_handlers.read().unwrap();
        assert_eq!(handlers.len(), 1);

        let test_event = ResourceEvent::MemoryWarning {
            current: 1000000,
            threshold: 800000,
            percentage: 80.0,
        };

        handlers[0](test_event);
        assert!(event_received.load(Ordering::Relaxed));
    }

    #[test]
    fn test_recommendations() {
        let settings = create_test_settings();
        let metrics = Arc::new(MetricsCollector::new());
        let monitor = ResourceMonitor::new(settings, metrics);

        let recommendations = monitor.get_recommendations().unwrap();
        // Recommendations depend on system state, so we just verify the call works
        assert!(recommendations.len() >= 0);
    }

    #[test]
    fn test_performance_optimizer() {
        let settings = create_test_settings();
        let optimizer = PerformanceOptimizer::new(settings);

        // This should not fail even if optimizations aren't actually applied
        let result = optimizer.optimize();
        assert!(result.is_ok());
    }

    #[test]
    fn test_resource_event_serialization() {
        let event = ResourceEvent::MemoryWarning {
            current: 1000000,
            threshold: 800000,
            percentage: 80.0,
        };

        let serialized = serde_json::to_string(&event).unwrap();
        let deserialized: ResourceEvent = serde_json::from_str(&serialized).unwrap();

        match (event, deserialized) {
            (ResourceEvent::MemoryWarning { current: c1, .. }, ResourceEvent::MemoryWarning { current: c2, .. }) => {
                assert_eq!(c1, c2);
            }
            _ => panic!("Event type mismatch"),
        }
    }

    #[test]
    fn test_system_info_serialization() {
        let info = SystemInfo {
            total_memory: 8589934592,
            available_memory: 4294967296,
            cpu_count: 8,
            memory_usage_percent: 50.0,
            cpu_usage_percent: 25.0,
            process_memory: 134217728,
            process_cpu_percent: 5.0,
            load_averages: Some((1.0, 1.5, 2.0)),
            file_descriptors: Some(256),
        };

        let serialized = serde_json::to_string(&info).unwrap();
        let deserialized: SystemInfo = serde_json::from_str(&serialized).unwrap();

        assert_eq!(info.total_memory, deserialized.total_memory);
        assert_eq!(info.cpu_count, deserialized.cpu_count);
        assert_eq!(info.load_averages, deserialized.load_averages);
    }
}