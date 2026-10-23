use std::fs;

/// Writes an exported file to the path the user picked in the native Save panel.
#[tauri::command]
fn write_file(path: String, contents: Vec<u8>) -> Result<(), String> {
    fs::write(&path, contents).map_err(|e| format!("Could not write {path}: {e}"))
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .invoke_handler(tauri::generate_handler![write_file])
        .run(tauri::generate_context!())
        .expect("error while running MolCanvas");
}
