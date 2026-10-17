//! edlib (vendored, rust/vendor/edlib), for placing a sequence within another allowing for
//! a set number of edits (its "HW" mode: free gaps at the target's ends).

use std::os::raw::{c_char, c_int, c_uchar};

#[repr(C)]
struct EdlibAlignConfig {
    k: c_int,
    mode: c_int,
    task: c_int,
    additional_equalities: *const c_char,
    additional_equalities_length: c_int,
}

#[repr(C)]
struct EdlibAlignResult {
    status: c_int,
    edit_distance: c_int,
    end_locations: *mut c_int,
    start_locations: *mut c_int,
    num_locations: c_int,
    alignment: *mut c_uchar,
    alignment_length: c_int,
    alphabet_length: c_int,
}

unsafe extern "C" {
    fn edlibAlign(
        query: *const c_char,
        query_length: c_int,
        target: *const c_char,
        target_length: c_int,
        config: EdlibAlignConfig,
    ) -> EdlibAlignResult;
    fn edlibFreeAlignResult(result: EdlibAlignResult);
}

const HW: c_int = 2;
const TASK_DISTANCE: c_int = 0;
const TASK_LOCATIONS: c_int = 1;

fn run(
    query: &[u8],
    target: &[u8],
    max_edits: usize,
    task: c_int,
) -> Option<(usize, Option<usize>)> {
    let config = EdlibAlignConfig {
        k: max_edits as c_int,
        mode: HW,
        task,
        additional_equalities: std::ptr::null(),
        additional_equalities_length: 0,
    };
    // SAFETY: the pointers and lengths come from live slices, and the result is freed below
    // after what's needed is copied out
    unsafe {
        let r = edlibAlign(
            query.as_ptr() as *const c_char,
            query.len() as c_int,
            target.as_ptr() as *const c_char,
            target.len() as c_int,
            config,
        );
        let out = if r.status != 0 || r.edit_distance < 0 {
            None
        } else if task == TASK_LOCATIONS {
            if r.num_locations < 1 || r.start_locations.is_null() {
                None
            } else {
                Some((r.edit_distance as usize, Some(*r.start_locations as usize)))
            }
        } else {
            Some((r.edit_distance as usize, None))
        };
        edlibFreeAlignResult(r);
        out
    }
}

/// (edits, start in target) of `query`'s first best placement within `target`, if within
/// max_edits
pub fn locate(query: &[u8], target: &[u8], max_edits: usize) -> Option<(usize, usize)> {
    run(query, target, max_edits, TASK_LOCATIONS).map(|(d, s)| (d, s.unwrap()))
}

/// edits of `query`'s best placement within `target`, if within max_edits
pub fn distance(query: &[u8], target: &[u8], max_edits: usize) -> Option<usize> {
    run(query, target, max_edits, TASK_DISTANCE).map(|(d, _)| d)
}
