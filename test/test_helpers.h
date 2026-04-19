#pragma once
#include <cstdlib>
#include <string>

// Returns the path to the LEEana working directory.
// Override with LEEANA_DIR env var; defaults to the known absolute path.
inline std::string leeana_dir() {
    const char* env = std::getenv("LEEANA_DIR");
    if (env && *env) return std::string(env);
    return "/home/xqian/wire-cell/wcp-uboone-bdt/LEEana";
}

inline std::string leeana_path(const std::string& rel) {
    return leeana_dir() + "/" + rel;
}
