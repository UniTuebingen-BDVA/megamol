/**
 * KeyframeManipulators.cpp
 *
 * Copyright (C) 2017 by VISUS (Universitaet Stuttgart).
 * Alle Rechte vorbehalten.
 */

#include "cinematic_gl/KeyframeManipulators.h"
#include "mmcore/utility/log/Log.h"
#include "vislib/math/Cuboid.h"
#include <glm/gtc/type_ptr.hpp>


using namespace megamol;
using namespace megamol::core;
using namespace megamol::cinematic_gl;


KeyframeManipulators::KeyframeManipulators(void)
        : toggleVisibleGroupParam(
              "manipulators::toggleVisibleGroup", "Toggle visibility of different manipulator groups.")
        , visibleGroupParam("manipulators::visibleGroup", "Select visible manipulator group.")
        , toggleOusideBboxParam(
              "manipulators::showOutsideBBox", "Show manipulators always outside of model bounding box.")
        , paramSlots()
        , visibleGroup(VisibleGroup::SELECTED_KEYFRAME_LOOKAT_AND_UP_VECTOR)
        , toggleOusideBbox(true)
        , manipulators()
        , selectors()
        , state() {

    this->paramSlots.clear();
    this->selectors.clear();
    this->manipulators.clear();

    this->toggleOusideBboxParam.SetParameter(
        new param::ButtonParam(core::view::Key::KEY_W, core::view::Modifier::SHIFT));
    this->paramSlots.emplace_back(&this->toggleOusideBboxParam);

    this->toggleVisibleGroupParam.SetParameter(
        new param::ButtonParam(core::view::Key::KEY_Q, core::view::Modifier::SHIFT));
    this->paramSlots.emplace_back(&this->toggleVisibleGroupParam);

    param::EnumParam* vmg = new param::EnumParam(this->visibleGroup);
    vmg->SetTypePair(VisibleGroup::SELECTED_KEYFRAME_AND_CTRLPOINT_POSITION, "Keyframe and Ctrl-Point Positions");
    vmg->SetTypePair(VisibleGroup::SELECTED_KEYFRAME_LOOKAT_AND_UP_VECTOR, "LookAt Vector and Up Vector");
    this->visibleGroupParam << vmg;
    this->paramSlots.emplace_back(&this->visibleGroupParam);
    vmg = nullptr;

    // Init state
    this->state.selected_keyframe = cinematic::Keyframe();
    this->state.viewport = glm::vec2(0.0f, 0.0f);
    this->state.mvp = glm::mat4();
    this->state.first_ctrl_point = glm::vec3(0.0f, 0.0f, 0.0f);
    this->state.last_ctrl_point = glm::vec3(0.0f, 0.0f, 0.0f);

    this->state.cam = view::Camera();
    auto intrinsics = core::view::Camera::PerspectiveParameters();
    intrinsics.fovy = 0.5f;
    intrinsics.aspect = 16.0f / 9.0f;
    intrinsics.near_plane = 0.01f;
    intrinsics.far_plane = 100.0f;
    /// intrinsics.image_plane_tile = ;
    this->state.cam.setPerspectiveProjection(intrinsics);

    this->state.bbox = vislib::math::Cuboid<float>(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    this->state.hit.reset();
    this->state.last_mouse = glm::vec2(0.0f, 0.0f);
    this->state.selected_index = 0;
    this->state.point_radius = 0.05f;
    this->state.line_width = 0.01f;
    this->state.line_length = 0.3f;

    // Init manipulators
    float lookat_length;
    Manipulator m;
    m.show = false;
    m.vector = glm::vec3(0.0f, 0.0f, 0.0f);

    m.variety = Manipulator::Variety::MANIPULATOR_FIRST_CTRLPOINT_POSITION;
    m.rigging = Manipulator::Rigging::X_DIRECTION;
    this->manipulators.push_back(m);
    m.rigging = Manipulator::Rigging::Y_DIRECTION;
    this->manipulators.push_back(m);
    m.rigging = Manipulator::Rigging::Z_DIRECTION;
    this->manipulators.push_back(m);

    m.variety = Manipulator::Variety::MANIPULATOR_LAST_CTRLPOINT_POSITION;
    m.rigging = Manipulator::Rigging::X_DIRECTION;
    this->manipulators.push_back(m);
    m.rigging = Manipulator::Rigging::Y_DIRECTION;
    this->manipulators.push_back(m);
    m.rigging = Manipulator::Rigging::Z_DIRECTION;
    this->manipulators.push_back(m);

    m.variety = Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_POSITION;
    m.rigging = Manipulator::Rigging::X_DIRECTION;
    this->manipulators.push_back(m);
    m.rigging = Manipulator::Rigging::Y_DIRECTION;
    this->manipulators.push_back(m);
    m.rigging = Manipulator::Rigging::Z_DIRECTION;
    this->manipulators.push_back(m);

    m.variety = Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_POSITION_USING_LOOKAT;
    m.rigging = Manipulator::Rigging::VECTOR_DIRECTION;
    this->manipulators.push_back(m);

    m.variety = Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_LOOKAT_VECTOR;
    m.rigging = Manipulator::Rigging::X_DIRECTION;
    this->manipulators.push_back(m);
    m.rigging = Manipulator::Rigging::Y_DIRECTION;
    this->manipulators.push_back(m);
    m.rigging = Manipulator::Rigging::Z_DIRECTION;
    this->manipulators.push_back(m);

    m.variety = Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_UP_VECTOR;
    m.rigging = Manipulator::Rigging::ROTATION;
    this->manipulators.push_back(m);
}


KeyframeManipulators::~KeyframeManipulators(void) {}


void KeyframeManipulators::UpdateExtents(vislib::math::Cuboid<float>& inout_bbox) {

    // First store copy of unmodified model bounding box
    this->state.bbox = inout_bbox;

    // Grow bounding box of model to manipulator positions
    glm::vec3 pos;
    for (auto& m : this->manipulators) {
        pos = this->getActualManipulatorPosition(m);
        inout_bbox.GrowToPoint(core_gl::utility::glm_to_vislib_point(pos));
    }
    for (auto& s : this->selectors) {
        pos = s.vector;
        pos = pos + glm::normalize(pos) * (this->state.point_radius * 2.0f);
        inout_bbox.GrowToPoint(core_gl::utility::glm_to_vislib_point(pos));
    }
}


bool KeyframeManipulators::UpdateRendering(const std::shared_ptr<std::vector<cinematic::Keyframe>> keyframes,
    cinematic::Keyframe selected_keyframe, glm::vec3 first_ctrl_pos, glm::vec3 last_ctrl_pos,
    core::view::Camera const& cam, glm::vec2 viewport_dim, glm::mat4 mvp) {

    // Update parameters
    if (this->visibleGroupParam.IsDirty()) {
        this->visibleGroup =
            static_cast<VisibleGroup>(this->visibleGroupParam.Param<core::param::EnumParam>()->Value());
        this->visibleGroupParam.ResetDirty();
    }

    if (this->toggleVisibleGroupParam.IsDirty()) {
        this->visibleGroup = static_cast<VisibleGroup>((this->visibleGroup + 1) % VisibleGroup::VISIBLEGROUP_COUNT);
        this->toggleVisibleGroupParam.ResetDirty();
    }

    if (this->toggleOusideBboxParam.IsDirty()) {
        this->toggleOusideBbox = !this->toggleOusideBbox;
        this->toggleOusideBboxParam.ResetDirty();
    }

    auto keyframe_count = keyframes->size();

    // Following stuff should not be updated during active manipulator dragging:
    if (this->state.hit == nullptr) {

        // Update current state
        this->state.viewport = viewport_dim;
        this->state.mvp = mvp;
        this->state.cam = cam;
        this->state.selected_keyframe = selected_keyframe;
        this->state.first_ctrl_point = first_ctrl_pos;
        this->state.last_ctrl_point = last_ctrl_pos;
        this->state.selected_index = -1;
        // Set primitive dimensions (empirical)
        float length = this->state.bbox.LongestEdge();
        this->state.line_length = length * 0.3f;
        this->state.line_width = length * 0.01f;
        this->state.point_radius = length * 0.05f;
        if (this->state.lookat_length <= 0.0f) {
            this->state.lookat_length = 2.0f * this->state.line_length;
        }

        // Update keyframe position selectors
        this->selectors.resize(keyframe_count);
        for (size_t i = 0; i < keyframe_count; ++i) {
            this->selectors[i].show = true;
            this->selectors[i].variety = Manipulator::Variety::SELECTOR_KEYFRAME_POSITION;
            this->selectors[i].rigging = Manipulator::Rigging::NONE;
            this->selectors[i].vector = keyframes->operator[](i).GetCamera().get<view::Camera::Pose>().position;
            if (keyframes->operator[](i) == this->state.selected_keyframe) {
                this->state.selected_index = i;
            }
        }
    }

    // Update selected keyframe manipulators
    view::Camera selected_camera = this->state.selected_keyframe.GetCamera();

    for (auto& m : this->manipulators) {

        // Set visibility of manipulators
        m.show = false;
        switch (this->visibleGroup) {
        case (VisibleGroup::SELECTED_KEYFRAME_AND_CTRLPOINT_POSITION): {
            switch (m.variety) {
            case (Manipulator::Variety::MANIPULATOR_FIRST_CTRLPOINT_POSITION): {
                if (this->state.selected_index == 0) {
                    m.show = true;
                }
            } break;
            case (Manipulator::Variety::MANIPULATOR_LAST_CTRLPOINT_POSITION): {
                if (this->state.selected_index == (keyframe_count - 1)) {
                    m.show = true;
                }
            } break;
            case (Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_POSITION): {
                m.show = true;
            } break;
            default:
                break;
            }
        } break;
        case (VisibleGroup::SELECTED_KEYFRAME_LOOKAT_AND_UP_VECTOR): {
            switch (m.variety) {
            case (Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_LOOKAT_VECTOR):
            case (Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_POSITION_USING_LOOKAT):
            case (Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_UP_VECTOR): {
                m.show = true;
            } break;
            default:
                break;
            }
        } break;
        default:
            break;
        }

        // Set direction vector of manipulators
        switch (m.rigging) {
        case (Manipulator::Rigging::X_DIRECTION): {
            m.vector = glm::vec3(1.0f, 0.0f, 0.0f);
        } break;
        case (Manipulator::Rigging::Y_DIRECTION): {
            m.vector = glm::vec3(0.0f, 1.0f, 0.0f);
        } break;
        case (Manipulator::Rigging::Z_DIRECTION): {
            m.vector = glm::vec3(0.0f, 0.0f, 1.0f);
        } break;
        case (Manipulator::Rigging::VECTOR_DIRECTION): {
            if (m.variety == Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_POSITION_USING_LOOKAT) {
                glm::vec3 view = selected_camera.get<view::Camera::Pose>().direction;
                m.vector = view * (-1.0f);
            }
        } break;
        case (Manipulator::Rigging::ROTATION): {
            if (m.variety == Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_UP_VECTOR) {
                m.vector = selected_camera.get<view::Camera::Pose>().up;
            }
        } break;
        default:
            break;
        }
        m.vector = glm::normalize(m.vector);
    }

    return true;
}


bool KeyframeManipulators::PushRendering(CinematicUtils& utils) {

    core::view::Camera global_camera(this->state.cam);
    auto global_cam_pose = global_camera.get<core::view::Camera::Pose>();
    glm::vec3 global_cam_position = global_cam_pose.position;
    glm::vec3 global_cam_view = global_cam_pose.direction;

    view::Camera selected_camera(this->state.selected_keyframe.GetCamera());
    auto selected_cam_pose = selected_camera.get<core::view::Camera::Pose>();
    glm::vec3 selected_cam_position = selected_cam_pose.position;
    glm::vec3 selected_cam_view = selected_cam_pose.direction;
    glm::vec3 selected_cam_up = selected_cam_pose.up;

    glm::vec3 manipulator_origin;
    glm::vec3 manipulator_position;
    glm::vec4 manipulator_color;

    // Push keyframe position selectors
    auto selector_count = this->selectors.size();
    for (size_t i = 0; i < selector_count; ++i) {
        manipulator_color = utils.Color(CinematicUtils::Colors::KEYFRAME);
        if (i == this->state.selected_index) {
            manipulator_color = utils.Color(CinematicUtils::Colors::KEYFRAME_SELECTED);
        }
        utils.PushPointPrimitive(this->selectors[i].vector, (2.0f * this->state.point_radius), global_cam_view,
            global_cam_position, manipulator_color);
        if (this->state.selected_index == i) {
            utils.PushLinePrimitive(this->selectors[i].vector, selected_cam_position, this->state.line_width,
                global_cam_view, global_cam_position, utils.Color(CinematicUtils::Colors::KEYFRAME_DRAGGED));
        }
    }

    if (this->state.selected_index >= 0) {
        // Push selected keyframe manipulators
        for (auto& m : this->manipulators) {
            if (m.show) {
                manipulator_color = this->getManipulatorColor(m, utils);
                manipulator_origin = this->getManipulatorOrigin(m);
                manipulator_position = this->getActualManipulatorPosition(m);
                switch (m.variety) {
                case (Manipulator::Variety::MANIPULATOR_FIRST_CTRLPOINT_POSITION):
                case (Manipulator::Variety::MANIPULATOR_LAST_CTRLPOINT_POSITION): {
                    if (m.rigging == Manipulator::Rigging::X_DIRECTION) { // Draw only once
                        utils.PushLinePrimitive(selected_cam_position, manipulator_origin, this->state.line_width,
                            global_cam_view, global_cam_position,
                            utils.Color(CinematicUtils::Colors::MANIPULATOR_CTRLPOINT));
                    }
                } break;
                case (Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_LOOKAT_VECTOR): {
                    if (m.rigging == Manipulator::Rigging::X_DIRECTION) { // Draw only once
                        utils.PushLinePrimitive(selected_cam_position, manipulator_origin, this->state.line_width,
                            global_cam_view, global_cam_position,
                            utils.Color(CinematicUtils::Colors::MANIPULATOR_VECTOR));
                    }
                } break;
                default:
                    break;
                }
                utils.PushPointPrimitive(manipulator_position, (2.0f * this->state.point_radius), global_cam_view,
                    global_cam_position, manipulator_color);
                utils.PushLinePrimitive(manipulator_origin, manipulator_position, this->state.line_width,
                    global_cam_view, global_cam_position, manipulator_color);
            }
        }
    } else {
        // Push intermediate keyframe marker
        selected_cam_view = glm::normalize(selected_cam_view) * this->state.line_length * 2.0f;
        selected_cam_up = glm::normalize(selected_cam_up) * this->state.line_length;
        manipulator_color = utils.Color(CinematicUtils::Colors::KEYFRAME);
        utils.PushPointPrimitive(selected_cam_position, (2.0f * this->state.point_radius) * (3.0f / 4.0f),
            global_cam_view, global_cam_position, manipulator_color);
        utils.PushLinePrimitive(selected_cam_position, selected_cam_position + selected_cam_view,
            this->state.line_width, global_cam_view, global_cam_position,
            utils.Color(CinematicUtils::Colors::MANIPULATOR_VECTOR));
        utils.PushLinePrimitive(selected_cam_position, selected_cam_position + selected_cam_up, this->state.line_width,
            global_cam_view, global_cam_position, utils.Color(CinematicUtils::Colors::MANIPULATOR_ROTATION));
    }

    return true;
}


int KeyframeManipulators::GetSelectedKeyframePositionIndex(float mouse_x, float mouse_y) {

    int index = -1;
    glm::vec2 mouse = glm::vec2(mouse_x, mouse_y);

    glm::vec3 cam_up = this->state.cam.get<core::view::Camera::Pose>().up;

    auto count = this->selectors.size();
    for (size_t i = 0; i < count; ++i) {
        if (this->selectors[i].show) {
            if (this->checkMousePointIntersection(this->selectors[i], mouse, cam_up)) {
                return i;
            }
        }
    }
    return index;
}


bool KeyframeManipulators::CheckForHitManipulator(float mouse_x, float mouse_y) {

    this->state.hit = nullptr;
    glm::vec2 mouse = glm::vec2(mouse_x, mouse_y);

    glm::vec3 cam_up = this->state.cam.get<core::view::Camera::Pose>().up;

    for (auto& m : this->manipulators) {
        if (m.show) {
            if (this->checkMousePointIntersection(m, mouse, cam_up)) {
                this->state.last_mouse = mouse;
                this->state.hit = std::make_shared<Manipulator>(m);
                return true;
            }
        }
    }
    return false;
}


bool KeyframeManipulators::ProcessHitManipulator(float mouse_x, float mouse_y) {

    if (this->state.hit == nullptr)
        return false;

    glm::vec3 manipulator_origin = this->getManipulatorOrigin((*this->state.hit));
    glm::vec3 manipulator_position = this->getActualManipulatorPosition((*this->state.hit));
    glm::vec3 world_vector = manipulator_position - manipulator_origin;
    float world_length = glm::length(world_vector);

    glm::vec2 screenspace_manipulator_origin = this->world2ScreenSpace(manipulator_origin);
    glm::vec2 screenspace_manipulator_position = this->world2ScreenSpace(manipulator_position);
    glm::vec2 screenspace_vector = screenspace_manipulator_position - screenspace_manipulator_origin;
    float screenspace_length = glm::length(screenspace_vector);

    glm::vec2 current_mouse = glm::vec2(mouse_x, mouse_y);
    glm::vec2 mouse_vector = current_mouse - this->state.last_mouse;
    screenspace_vector = glm::normalize(screenspace_vector);
    float diff_screenspace_length = glm::dot(screenspace_vector, mouse_vector);

    float diff_world_length = world_length / screenspace_length * diff_screenspace_length;
    glm::vec3 diff_world_vector = this->state.hit->vector * diff_world_length;

    auto selected_camera = this->state.selected_keyframe.GetCamera();
    auto selected_camera_pose = selected_camera.get<core::view::Camera::Pose>();
    glm::vec3 camera_position = selected_camera_pose.position;

    if (this->state.hit->variety == Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_LOOKAT_VECTOR) {
        // Handle position of look at manipulator
        glm::vec3 new_view = (manipulator_origin - camera_position) + diff_world_vector;
        selected_camera_pose.direction = new_view;

        this->state.lookat_length = glm::length(new_view);
    } else if (this->state.hit->variety == Manipulator::Variety::MANIPULATOR_FIRST_CTRLPOINT_POSITION) {
        // Handle position of first control point manipulator
        this->state.first_ctrl_point = manipulator_origin + diff_world_vector;
    } else if (this->state.hit->variety == Manipulator::Variety::MANIPULATOR_LAST_CTRLPOINT_POSITION) {
        // Handle position of last control point manipulator
        this->state.last_ctrl_point = manipulator_origin + diff_world_vector;
    } else if ((this->state.hit->variety == Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_UP_VECTOR) &&
               (this->state.hit->rigging == Manipulator::Rigging::ROTATION)) {
        // Handle rotation of up vector
        glm::vec3 view = selected_camera_pose.direction;
        glm::vec3 right = glm::cross(selected_camera_pose.direction, selected_camera_pose.up);

        // Using right vector as delta vector for manipulating up vector rotating around view vector
        glm::vec3 right_position = manipulator_origin + glm::normalize(right) * world_length;
        glm::vec2 screenspace_right_position = this->world2ScreenSpace(right_position);
        screenspace_vector = screenspace_right_position - screenspace_manipulator_origin;
        screenspace_length = glm::length(screenspace_vector);
        screenspace_vector = glm::normalize(screenspace_vector);
        diff_screenspace_length = glm::dot(screenspace_vector, mouse_vector);
        diff_world_length = world_length / screenspace_length * diff_screenspace_length;
        diff_world_vector = right * diff_world_length;

        glm::vec3 new_up = world_vector + diff_world_vector;
        selected_camera_pose.up = new_up;

        this->state.hit->vector = glm::normalize(new_up);
    } else {
        glm::vec3 new_camera_position =
            camera_position + glm::vec3(diff_world_vector.x, diff_world_vector.y, diff_world_vector.z);
        selected_camera_pose.position = new_camera_position;
    }

    // Apply changed camera state to selected keyframe
    selected_camera_pose.direction = glm::normalize(selected_camera_pose.direction);
    selected_camera_pose.up = glm::normalize(selected_camera_pose.up);
    selected_camera_pose.right = glm::normalize(selected_camera_pose.right);
    selected_camera.setPose(selected_camera_pose);
    this->state.selected_keyframe.SetCameraState(selected_camera);

    this->state.last_mouse = current_mouse;

    return true;
}


glm::vec3 KeyframeManipulators::getManipulatorOrigin(Manipulator& manipulator) {

    glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f);
    switch (manipulator.variety) {
    case (Manipulator::Variety::SELECTOR_KEYFRAME_POSITION): {
        position = manipulator.vector;
    } break;
    case (Manipulator::Variety::MANIPULATOR_FIRST_CTRLPOINT_POSITION): {
        position = this->state.first_ctrl_point;
    } break;
    case (Manipulator::Variety::MANIPULATOR_LAST_CTRLPOINT_POSITION): {
        position = this->state.last_ctrl_point;
    } break;
    case (Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_POSITION): {
        position = this->state.selected_keyframe.GetCamera().get<view::Camera::Pose>().position;
    } break;
    case (Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_POSITION_USING_LOOKAT): {
        position = this->state.selected_keyframe.GetCamera().get<view::Camera::Pose>().position;
    } break;
    case (Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_LOOKAT_VECTOR): {
        position = this->state.selected_keyframe.GetCamera().get<view::Camera::Pose>().position;
        view::Camera cam(this->state.selected_keyframe.GetCamera());
        glm::vec3 snap_view = cam.get<view::Camera::Pose>().direction;
        snap_view = glm::normalize(snap_view);
        position = position + glm::vec3(snap_view.x, snap_view.y, snap_view.z) * this->state.lookat_length;
    } break;
    case (Manipulator::Variety::MANIPULATOR_SELECTED_KEYFRAME_UP_VECTOR): {
        position = this->state.selected_keyframe.GetCamera().get<view::Camera::Pose>().position;
    } break;
    default:
        break;
    }
    return position;
}


glm::vec3 KeyframeManipulators::getActualManipulatorPosition(Manipulator& manipulator) {

    glm::vec3 manipulator_origin = this->getManipulatorOrigin(manipulator);
    glm::vec3 manipulator_position = manipulator_origin + (manipulator.vector * this->state.line_length);

    if (this->toggleOusideBbox && (this->state.line_length > 0.0f) && (glm::length(manipulator.vector) > 0.0f) &&
        (!this->state.bbox.IsEmpty())) {
        float offset = 0.0f;
        while (this->state.bbox.Contains(core_gl::utility::glm_to_vislib_point(manipulator_position))) {
            offset += this->state.line_length / 4.0f;
            manipulator_position += manipulator.vector * offset;
        }
    }

    return manipulator_position;
}


glm::vec4 KeyframeManipulators::getManipulatorColor(Manipulator& manipulator, CinematicUtils& utils) {

    glm::vec4 color = glm::vec4(0.0f, 0.0f, 0.0f, 0.0f);
    switch (manipulator.rigging) {
    case (Manipulator::Rigging::X_DIRECTION): {
        color = utils.Color(CinematicUtils::Colors::MANIPULATOR_X);
    } break;
    case (Manipulator::Rigging::Y_DIRECTION): {
        color = utils.Color(CinematicUtils::Colors::MANIPULATOR_Y);
    } break;
    case (Manipulator::Rigging::Z_DIRECTION): {
        color = utils.Color(CinematicUtils::Colors::MANIPULATOR_Z);
    } break;
    case (Manipulator::Rigging::VECTOR_DIRECTION): {
        color = utils.Color(CinematicUtils::Colors::MANIPULATOR_VECTOR);
    } break;
    case (Manipulator::Rigging::ROTATION): {
        color = utils.Color(CinematicUtils::Colors::MANIPULATOR_ROTATION);
    } break;
    default:
        break;
    }
    return color;
}


bool KeyframeManipulators::checkMousePointIntersection(Manipulator& manipulator, glm::vec2 mouse, glm::vec3 cam_up) {

    glm::vec3 manipulator_position;
    if (manipulator.variety == Manipulator::Variety::SELECTOR_KEYFRAME_POSITION) {
        manipulator_position = this->getManipulatorOrigin(manipulator);
    } else {
        manipulator_position = this->getActualManipulatorPosition(manipulator);
    }
    glm::vec3 manipulator_extent = manipulator_position + glm::normalize(cam_up) * this->state.point_radius;

    glm::vec2 position_screenspace = this->world2ScreenSpace(manipulator_position);
    glm::vec2 extent_screenspace = this->world2ScreenSpace(manipulator_extent);
    float radius_screenspace = glm::length(position_screenspace - extent_screenspace);

    return (glm::length(position_screenspace - mouse) <= radius_screenspace);
}


glm::vec2 KeyframeManipulators::world2ScreenSpace(glm::vec3 vec) {

    glm::vec4 world = {vec.x, vec.y, vec.z, 1.0f};
    world = this->state.mvp * world;
    world = world / world.w;
    glm::vec2 screen;
    screen.x = (world.x + 1.0f) / 2.0f * this->state.viewport.x;
    screen.y = glm::abs(world.y - 1.0f) / 2.0f * this->state.viewport.y; // (flipped y-axis)

    return screen;
}
