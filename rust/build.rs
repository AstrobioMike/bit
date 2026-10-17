// compiles the vendored edlib (C++), used for the edit-distance alignments in find-overrep-seqs
fn main() {
    cc::Build::new()
        .cpp(true)
        .file("vendor/edlib/edlib.cpp")
        .include("vendor/edlib")
        .flag_if_supported("-std=c++11")
        .opt_level(3)
        .warnings(false)
        .compile("edlib");
    println!("cargo:rerun-if-changed=vendor/edlib/edlib.cpp");
    println!("cargo:rerun-if-changed=vendor/edlib/edlib.h");
}
