// Test script to understand sufr API
use sufr::sa::SuffixArray;

fn main() {
    let text = b"banana";
    
    // Try to create suffix array with sufr
    match SuffixArray::new(text) {
        Ok(sa) => {
            println!("Suffix array created successfully");
            println!("Length: {}", sa.len());
            
            // Try to get the actual suffix array data
            for i in 0..sa.len() {
                println!("sa[{}] = {}", i, sa.sa()[i]);
            }
        }
        Err(e) => {
            println!("Error creating suffix array: {:?}", e);
        }
    }
}