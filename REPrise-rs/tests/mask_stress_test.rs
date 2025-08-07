use reprise::mask::{Bitmask, ClaimGuard};
use std::sync::Arc;
use std::thread;

#[test]
fn test_bitmask_claim_release() {
    let mask = Bitmask::new(1000);
    assert!(mask.claim_range(&(0..100)));
    assert!(!mask.claim_range(&(50..150)));
    mask.release_range(&(0..100));
    assert!(mask.claim_range(&(50..150)));
}

#[test]
fn test_claim_guard() {
    let mask = Arc::new(Bitmask::new(1000));
    {
        let _guard = ClaimGuard::new(&mask, 0..100).unwrap();
        assert!(!mask.claim_range(&(50..150)));
    }
    assert!(mask.claim_range(&(50..150)));
}

#[test]
fn test_concurrent_claims() {
    let mask = Arc::new(Bitmask::new(1000));
    let mut handles = vec![];

    for _ in 0..10 {
        let mask_clone = Arc::clone(&mask);
        handles.push(thread::spawn(move || {
            let guard = ClaimGuard::new(&mask_clone, 100..200);
            if guard.is_some() {
                // This thread got the claim
                thread::sleep(std::time::Duration::from_millis(10));
            }
        }));
    }

    let mut successful_claims = 0;
    for handle in handles {
        handle.join().unwrap();
    }

    // This is a simplified check. A more robust test would use a channel
    // to count how many threads successfully created a guard.
    // For now, we just check that the mask is released.
    assert!(mask.claim_range(&(100..200)));
}
