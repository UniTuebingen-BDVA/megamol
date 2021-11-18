/*
 * ImagePresentationEntryPoints.h
 *
 * Copyright (C) 2021 by VISUS (Universitaet Stuttgart).
 * Alle Rechte vorbehalten.
 */

#pragma once

#include "EntryPoint.h"
#include "FrontendResource.h"

#include <functional>
#include <vector>
#include <any>

namespace megamol {
namespace frontend_resources {

using EntryPointRenderFunctions = std::tuple<
    // ptr to entry point object
    void*,
    // rendering execution function
    EntryPointExecutionCallback,
    // get requested resources function
    std::function<std::vector<std::string>()>>;

struct ImagePresentationEntryPoints {
    std::function<bool(std::string const&, EntryPointRenderFunctions const&)> add_entry_point;
    std::function<bool(std::string const&)> remove_entry_point;
    std::function<bool(std::string const&, std::string const&)> rename_entry_point;
    std::function<void()> clear_entry_points;

    enum class SubscriptionEvent {
        Add, Remove, Rename, Clear
    };

    using SubscriberFunction = std::function<
        void(
            frontend_resources::ImagePresentationEntryPoints::SubscriptionEvent const&,
            std::vector<std::any> const& /*args*/)
    >;

    std::function<void( SubscriberFunction const& )> subscribe_to_entry_point_changes;
    std::function<optional<EntryPoint>(std::string const&)> get_entry_point;
};

} /* end namespace frontend_resources */
} /* end namespace megamol */
