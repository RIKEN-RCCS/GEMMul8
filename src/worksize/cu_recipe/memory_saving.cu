#include "../../config/config.hpp"

#include <atomic>
#include <memory>
#include <mutex>
#include <unordered_map>

namespace gemmul8::config {

namespace {

struct Config {
    std::atomic<bool> memory_saving{DEFAULT_MEMORY_SAVING};
    std::atomic<size_t> max_worksize{DEFAULT_MAX_WORKSIZE};
    std::atomic<int> block_size_trsm{DEFAULT_BLOCK_SIZE_TRSM};
};

using ConfigPtr = std::shared_ptr<Config>;

struct BlasRegistryTag {};
struct LtRegistryTag {};

template <typename Tag, typename Handle>
struct ConfigRegistry {
    std::mutex mutex;
    std::unordered_map<Handle, ConfigPtr> map;
};

template <typename Tag, typename Handle>
ConfigRegistry<Tag, Handle> &config_registry() noexcept {
    static ConfigRegistry<Tag, Handle> registry;
    return registry;
}

template <typename Tag, typename Handle>
ConfigPtr get_or_create_config_ptr(Handle handle) noexcept {
    if (handle == nullptr) return {};

    auto &registry = config_registry<Tag, Handle>();
    try {
        std::lock_guard<std::mutex> lock(registry.mutex);
        auto &ptr = registry.map[handle];
        if (!ptr) ptr = std::make_shared<Config>();
        return ptr;
    } catch (...) {
        return {};
    }
}

template <typename Tag, typename Handle>
ConfigPtr find_config_ptr(Handle handle) noexcept {
    if (handle == nullptr) return {};

    auto &registry = config_registry<Tag, Handle>();
    try {
        std::lock_guard<std::mutex> lock(registry.mutex);
        const auto it = registry.map.find(handle);
        return (it == registry.map.end()) ? ConfigPtr{} : it->second;
    } catch (...) {
        return {};
    }
}

template <typename Tag, typename Handle>
void set_config_ptr(Handle handle, const ConfigPtr &ptr) noexcept {
    if (handle == nullptr || !ptr) return;

    auto &registry = config_registry<Tag, Handle>();
    try {
        std::lock_guard<std::mutex> lock(registry.mutex);
        registry.map[handle] = ptr;
    } catch (...) {
    }
}

template <typename Tag, typename Handle>
void clear_config_t(Handle handle) noexcept {
    if (handle == nullptr) return;

    auto &registry = config_registry<Tag, Handle>();
    try {
        std::lock_guard<std::mutex> lock(registry.mutex);
        registry.map.erase(handle);
    } catch (...) {
    }
}

ConfigSnapshot snapshot_config(const ConfigPtr &ptr) noexcept {
    if (!ptr) return {};

    ConfigSnapshot out;
    out.memory_saving   = ptr->memory_saving.load(std::memory_order_relaxed);
    out.max_worksize    = ptr->max_worksize.load(std::memory_order_relaxed);
    out.block_size_trsm = ptr->block_size_trsm.load(std::memory_order_relaxed);
    return out;
}

template <typename Tag, typename Handle, typename F>
void update_config(Handle handle, F &&update) noexcept {
    if (handle == nullptr) return;

    auto &registry = config_registry<Tag, Handle>();
    try {
        std::lock_guard<std::mutex> lock(registry.mutex);
        auto &ptr = registry.map[handle];
        if (!ptr) ptr = std::make_shared<Config>();
        update(*ptr);
    } catch (...) {
    }
}

} // namespace

ConfigSnapshot get_config(cublasHandle_t handle) noexcept {
    return snapshot_config(find_config_ptr<BlasRegistryTag>(handle));
}

ConfigSnapshot get_configLt(cublasLtHandle_t handle) noexcept {
    return snapshot_config(find_config_ptr<LtRegistryTag>(handle));
}

void set_memory_saving_impl(cublasHandle_t handle, bool enable) noexcept {
    update_config<BlasRegistryTag>(handle, [&](Config &config) {
        config.memory_saving.store(enable, std::memory_order_relaxed);
    });
}

void set_memory_savingLt_impl(cublasLtHandle_t handle, bool enable) noexcept {
    update_config<LtRegistryTag>(handle, [&](Config &config) {
        config.memory_saving.store(enable, std::memory_order_relaxed);
    });
}

void set_max_worksize_impl(cublasHandle_t handle, size_t bytes) noexcept {
    update_config<BlasRegistryTag>(handle, [&](Config &config) {
        config.max_worksize.store(bytes, std::memory_order_relaxed);
    });
}

void set_max_worksizeLt_impl(cublasLtHandle_t handle, size_t bytes) noexcept {
    update_config<LtRegistryTag>(handle, [&](Config &config) {
        config.max_worksize.store(bytes, std::memory_order_relaxed);
    });
}

void set_block_size_trsm_impl(cublasHandle_t handle, int nB) noexcept {
    update_config<BlasRegistryTag>(handle, [&](Config &config) {
        config.block_size_trsm.store(nB, std::memory_order_relaxed);
    });
}

void set_block_size_trsmLt_impl(cublasLtHandle_t handle, int nB) noexcept {
    update_config<LtRegistryTag>(handle, [&](Config &config) {
        config.block_size_trsm.store(nB, std::memory_order_relaxed);
    });
}

void clear_config_impl(cublasHandle_t handle) noexcept {
    clear_config_t<BlasRegistryTag>(handle);
}

void clear_configLt_impl(cublasLtHandle_t handle) noexcept {
    clear_config_t<LtRegistryTag>(handle);
}

void bind_config(cublasHandle_t parent, cublasLtHandle_t child) noexcept {
    const auto ptr = get_or_create_config_ptr<BlasRegistryTag>(parent);
    set_config_ptr<LtRegistryTag>(child, ptr);
}

void unbind_configLt(cublasLtHandle_t child) noexcept {
    clear_config_t<LtRegistryTag>(child);
}

} // namespace gemmul8::config

namespace gemmul8 {

void set_memory_saving(cublasHandle_t handle, bool enable) noexcept {
    config::set_memory_saving_impl(handle, enable);
}

void set_memory_savingLt(cublasLtHandle_t handle, bool enable) noexcept {
    config::set_memory_savingLt_impl(handle, enable);
}

bool get_memory_saving(cublasHandle_t handle) noexcept {
    return config::get_config(handle).memory_saving;
}

bool get_memory_savingLt(cublasLtHandle_t handle) noexcept {
    return config::get_configLt(handle).memory_saving;
}

void set_max_worksize(cublasHandle_t handle, size_t bytes) noexcept {
    config::set_max_worksize_impl(handle, bytes);
}

void set_max_worksizeLt(cublasLtHandle_t handle, size_t bytes) noexcept {
    config::set_max_worksizeLt_impl(handle, bytes);
}

size_t get_max_worksize(cublasHandle_t handle) noexcept {
    return config::get_config(handle).max_worksize;
}

size_t get_max_worksizeLt(cublasLtHandle_t handle) noexcept {
    return config::get_configLt(handle).max_worksize;
}

void clear_config(cublasHandle_t handle) noexcept {
    config::clear_config_impl(handle);
}

void clear_configLt(cublasLtHandle_t handle) noexcept {
    config::clear_configLt_impl(handle);
}

} // namespace gemmul8
