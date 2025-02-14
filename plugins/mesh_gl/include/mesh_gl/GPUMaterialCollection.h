/*
 * GPURenderTaskCollection.h
 *
 * Copyright (C) 2019 by Universitaet Stuttgart (VISUS).
 * All rights reserved.
 */

#ifndef GPU_MATERIAL_COLLECTION_H_INCLUDED
#define GPU_MATERIAL_COLLECTION_H_INCLUDED

#include <memory>
#include <unordered_map>
#include <variant>
#include <vector>

#include <glowl/glowl.h>

#include "mmcore/CoreInstance.h"

namespace megamol {
namespace mesh_gl {

typedef glowl::GLSLProgram Shader;

class GPUMaterialCollection {
public:
    using TexturePtrType = std::variant<std::shared_ptr<glowl::Texture>, std::shared_ptr<glowl::Texture2D>,
        std::shared_ptr<glowl::Texture2DArray>, std::shared_ptr<glowl::Texture3D>,
        std::shared_ptr<glowl::TextureCubemapArray>>;

    struct Material {
        std::shared_ptr<Shader> shader_program;
        std::vector<std::shared_ptr<glowl::Texture>> textures;
    };

    void addMaterial(megamol::core::CoreInstance* mm_core_inst, std::string const& identifier,
        std::vector<std::filesystem::path> const& shader_filepaths,
        std::vector<std::shared_ptr<glowl::Texture>> const& textures = {});

    void addMaterial(std::string const& identifier, std::shared_ptr<Shader> const& shader,
        std::vector<std::shared_ptr<glowl::Texture>> const& textures = {});

    void addMaterial(std::string const& identifier, Material const& material);

    void updateMaterialTexture(
        std::string const& identifier, size_t tex_idx, std::shared_ptr<glowl::Texture> const& texture);

    void deleteMaterial(std::string const& identifier);

    void clear();

    Material const& getMaterial(std::string const& identifier);

    inline std::unordered_map<std::string, GPUMaterialCollection::Material> const& getMaterials() {
        return m_materials;
    }

private:
    std::unordered_map<std::string, Material> m_materials;
};

} // namespace mesh_gl
} // namespace megamol

#endif // !GPU_MATERIAL_COLLECTION_H_INCLUDED
