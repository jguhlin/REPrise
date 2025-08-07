// Debug test to check paths
#[test]
fn debug_working_directory() {
    let cwd = std::env::current_dir().unwrap();
    println!("Current working directory: {}", cwd.display());
    
    let cpp_path = cwd.parent().unwrap().join("REPrise");
    println!("C++ executable path: {}", cpp_path.display());
    println!("C++ exists: {}", cpp_path.exists());
    
    let ecoli_path = cwd.parent().unwrap().join("data").join("ecoli.fasta");
    println!("E. coli file path: {}", ecoli_path.display());
    println!("E. coli exists: {}", ecoli_path.exists());
}