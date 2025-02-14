/*
 * thecam/translate_manipulator.h
 *
 * Copyright (C) 2016 TheLib Team (http://www.thelib.org/license)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * - Neither the name of TheLib, TheLib Team, nor the names of its
 *   contributors may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THELIB TEAM AS IS AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL THELIB TEAM BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/*
 * megamol::core::thecam::translate_manipulator<T>::translate_manipulator
 */
template<class T>
megamol::core::thecam::translate_manipulator<T>::translate_manipulator(const world_type stepSize)
        : stepSize(stepSize) {}


/*
 * megamol::core::thecam::translate_manipulator<T>::~translate_manipulator
 */
template<class T>
megamol::core::thecam::translate_manipulator<T>::~translate_manipulator(void) {}


/*
 *  megamol::core::thecam::translate_manipulator<T>::move_forward
 */
template<class T>
void megamol::core::thecam::translate_manipulator<T>::move_forward(const world_type dist) {
    if (this->enabled()) {
        auto cam = this->camera();
        auto cam_pose = cam->template get<view::Camera::Pose>();
        cam_pose.position += dist * cam_pose.direction;
        cam->setPose(cam_pose);
    }
}


/*
 * megamol::core::thecam::translate_manipulator<T>::move_horizontally
 */
template<class T>
void megamol::core::thecam::translate_manipulator<T>::move_horizontally(const world_type dist) {
    if (this->enabled()) {
        auto cam = this->camera();
        auto cam_pose = cam->template get<view::Camera::Pose>();
        auto cam_right = glm::cross(cam_pose.direction, cam_pose.up);
        cam_pose.position += dist * cam_right;
        cam->setPose(cam_pose);
    }
}


/*
 * megamol::core::thecam::translate_manipulator<T>::move_vertically
 */
template<class T>
void megamol::core::thecam::translate_manipulator<T>::move_vertically(const world_type dist) {
    if (this->enabled()) {
        auto cam = this->camera();
        auto cam_pose = cam->template get<view::Camera::Pose>();
        cam_pose.position += dist * cam_pose.up;
        cam->setPose(cam_pose);
    }
}
