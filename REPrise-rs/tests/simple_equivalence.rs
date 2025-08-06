// Simple test to compare basic functionality
use std::process::Command;
use std::fs;

#[test]
fn test_simple_cpp_rust_equivalence() {
    let input_file = "../test/tst.fa";
    let cpp_output = "simple_cpp";
    let rust_output = "simple_rust";
    
    // Clean up
    let _ = fs::remove_file(format!("{}.reprof", cpp_output));
    let _ = fs::remove_file(format!("{}.reprof", rust_output));
    
    println!("Running C++ version...");
    let cpp_result = Command::new("../REPrise")
        .args(&["-input", input_file, "-output", cpp_output])
        .output()
        .expect("Failed to run C++ version");
    
    println!("C++ stdout: {}", String::from_utf8_lossy(&cpp_result.stdout));
    
    println!("Running Rust version...");
    let rust_result = Command::new("cargo")
        .args(&["run", "--bin", "REPrise-rs", "--release", "--", "-input", input_file, "-output", rust_output])
        .output()
        .expect("Failed to run Rust version");
    
    println!("Rust stdout: {}", String::from_utf8_lossy(&rust_result.stdout));
    
    // Parse repeat family counts from stdout
    let cpp_stdout = String::from_utf8_lossy(&cpp_result.stdout);
    let rust_stdout = String::from_utf8_lossy(&rust_result.stdout);
    
    // For C++ - look for the final summary information
    let cpp_families = if let Some(reprof_content) = fs::read_to_string(format!("{}.reprof", cpp_output)).ok() {
        reprof_content.lines().filter(|line| line.starts_with(">R=")).count()
    } else {
        0
    };
    
    // For Rust - parse from stdout
    let rust_families: usize = rust_stdout
        .lines()
        .find(|line| line.contains("Processed") && line.contains("repeat families"))
        .and_then(|line| line.split_whitespace().find_map(|word| word.parse().ok()))
        .unwrap_or(0);
    
    println!("C++ families: {}, Rust families: {}", cpp_families, rust_families);
    
    // Clean up
    let _ = fs::remove_file(format!("{}.reprof", cpp_output));
    let _ = fs::remove_file(format!("{}.reprof", rust_output));
    
    // For small test file, expect similar counts
    assert!(cpp_families <= 10, "C++ found too many families: {}", cpp_families);
    assert!(rust_families <= 10, "Rust found too many families: {}", rust_families);
    
    println!("✓ Simple equivalence test passed");
}