// Test sufr API
fn main() {
    println!("Testing sufr crate...");
    
    // Try to use the sufr functions
    let text = b"banana";
    
    // Check what functions are available
    println!("Text: {:?}", std::str::from_utf8(text).unwrap());
    
    // Try to call sufr::create or similar
    // sufr::create(...);
}