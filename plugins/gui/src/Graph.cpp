/*
 * Graph.cpp
 *
 * Copyright (C) 2019 by Universitaet Stuttgart (VIS).
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "Graph.h"


using namespace megamol;
using namespace megamol::gui;
using vislib::sys::Log;


megamol::gui::Graph::Graph(void)
    : modules_graph()
    , calls_graph()
    , modules_list()
    , calls_list() {

}


megamol::gui::Graph::~Graph(void) {

}


bool megamol::gui::Graph::AddModule(const std::string& module_class_name) {

    Graph::Module mod;
    mod.class_name.clear();
    for (auto& m : this->modules_list) {
        if (module_class_name == m.class_name) {
            mod = m;
        }
    }
    if (mod.class_name.empty()) {
        vislib::sys::Log::DefaultLog.WriteError("Unable to find module for class name: %s [%s, %s, line %d]\n", module_class_name.c_str(), __FILE__, __FUNCTION__, __LINE__);
        return false;
    }

///TODO Replace dummy names
    mod.name = mod.class_name + "XXX";
    mod.full_name = "XXX::" + mod.name;
    mod.instance = "inst";

///TODO Adjuist size depending on size of content (module name, slot name)
    mod.gui.position = ImVec2(10.0f, 10.0f);
    mod.gui.size = ImVec2(200.0f, 100.0f);

    this->modules_graph.emplace_back(mod);

    return true;
}


bool megamol::gui::Graph::AddCall(const std::string& call_class_name) {


    return true;
}


bool megamol::gui::Graph::UpdateAvailableModulesCallsOnce(const megamol::core::CoreInstance* core_instance) {

    bool retval = true;
    if (core_instance == nullptr) {
        vislib::sys::Log::DefaultLog.WriteError("Pointer to Core Instance is nullptr. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        return false;
    }

    // CALLS ------------------------------------------------------------------
    /// ! Get calls before getting modules for having calls in place for setting compatible call indices of slots!
    if (this->calls_list.empty()) {

        // Get core calls
        std::string plugin_name = "Core";
        for (auto& c_desc : core_instance->GetCallDescriptionManager()) {
            Graph::Call call;
            call.plugin_name = plugin_name;
            this->read_call_data(call, c_desc);
            calls_list.emplace_back(call);
        }
        //Get plugin calls
        const std::vector<core::utility::plugins::AbstractPluginInstance::ptr_type>& plugins = core_instance->Plugins().GetPlugins();
        for (core::utility::plugins::AbstractPluginInstance::ptr_type plugin : plugins) {
            plugin_name = plugin->GetAssemblyName();
            for (auto& c_desc : plugin->GetCallDescriptionManager()) {
                Graph::Call call;
                call.plugin_name = plugin_name;
                this->read_call_data(call, c_desc);
                calls_list.emplace_back(call);
            }
        }

        vislib::sys::Log::DefaultLog.WriteInfo("Reading available calls ... DONE. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
    }

    // MODULES ----------------------------------------------------------------
    if (this->modules_list.empty()) {
///TODO Paralllalize using threads or MPI?
        // Allocate memory for each module and pass per reference to thread;
        /*
        auto module_count = core_instance->GetModuleDescriptionManager().Count();
        const std::vector<core::utility::plugins::AbstractPluginInstance::ptr_type>& plugins = core_instance->Plugins().GetPlugins();
        for (core::utility::plugins::AbstractPluginInstance::ptr_type plugin : plugins) {
            module_count += plugin->GetModuleDescriptionManager().Count();
        }
        this->modules_list.resize(module_count);
        */

        // Get core modules
        std::string plugin_name = "Core";
        for (auto& m_desc : core_instance->GetModuleDescriptionManager()) {
            Graph::Module mod;
            mod.plugin_name = std::string("core");
            this->read_module_data(mod, m_desc, core_instance);
            this->modules_list.emplace_back(mod);
        }
        // Get plugin modules
        const std::vector<core::utility::plugins::AbstractPluginInstance::ptr_type>& plugins = core_instance->Plugins().GetPlugins();
        for (core::utility::plugins::AbstractPluginInstance::ptr_type plugin : plugins) {
            plugin_name = plugin->GetAssemblyName();
            for (auto& m_desc : plugin->GetModuleDescriptionManager()) {
                Graph::Module mod;
                mod.plugin_name = plugin_name;
                this->read_module_data(mod, m_desc, core_instance);
                this->modules_list.emplace_back(mod);
            }
        }
        // Sorting module names alphabetically ascending
        std::sort(this->modules_list.begin(), this->modules_list.end(), [](Graph::Module mod1, Graph::Module mod2) {
            std::vector<std::string> v;
            v.clear();
            v.emplace_back(mod1.class_name);
            v.emplace_back(mod2.class_name);
            std::sort(v.begin(), v.end());
            return (v.front() != mod2.class_name);
        });

        vislib::sys::Log::DefaultLog.WriteInfo("Reading available modules ... DONE. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
    }

    return retval;
}


bool megamol::gui::Graph::read_module_data(Graph::Module& mod, const std::shared_ptr<const megamol::core::factories::ModuleDescription> mod_desc, const megamol::core::CoreInstance* core_instance) {

    mod.class_name = std::string(mod_desc->ClassName());
    mod.description = std::string(mod_desc->Description());
    /// mod.plugin_name is not available here (set above). 
    mod.param_slots.clear();
    mod.callee_slots.clear();
    mod.caller_slots.clear();

    // SLOTS ------------------------------------------------------------------
    /// (Following code is adapted from megamol::core::job::job::PluginsStateFileGeneratorJob.cpp)

    megamol::core::Module::ptr_type new_mod(mod_desc->CreateModule(nullptr));
    if (new_mod == nullptr) {
        vislib::sys::Log::DefaultLog.WriteError("Unable to create module: %s. [%s, %s, line %d]\n", mod_desc->ClassName(), __FILE__, __FUNCTION__, __LINE__);
        return false;
    }
    megamol::core::RootModuleNamespace::ptr_type root_mod_ns = std::make_shared<megamol::core::RootModuleNamespace>();
    root_mod_ns->AddChild(new_mod);

    std::shared_ptr<const core::view::AbstractView> viewptr = std::dynamic_pointer_cast<const core::view::AbstractView>(new_mod);
    mod.is_view = (viewptr != nullptr);

    std::vector<std::shared_ptr<core::param::ParamSlot> > paramSlots;
    std::vector<std::shared_ptr<core::CallerSlot> > callerSlots;
    std::vector<std::shared_ptr<core::CalleeSlot> > calleeSlots;

    core::Module::child_list_type::iterator ano_end = new_mod->ChildList_End();
    for (core::Module::child_list_type::iterator ano_i = new_mod->ChildList_Begin(); ano_i != ano_end; ++ano_i) {
        std::shared_ptr<core::param::ParamSlot> p_ptr = std::dynamic_pointer_cast<core::param::ParamSlot>(*ano_i);
        if (p_ptr != nullptr) paramSlots.push_back(p_ptr);
        std::shared_ptr<core::CallerSlot> cr_ptr = std::dynamic_pointer_cast<core::CallerSlot>(*ano_i);
        if (cr_ptr != nullptr) callerSlots.push_back(cr_ptr);
        std::shared_ptr<core::CalleeSlot> ce_ptr = std::dynamic_pointer_cast<core::CalleeSlot>(*ano_i);
        if (ce_ptr != nullptr) calleeSlots.push_back(ce_ptr);
    }

    // Param Slots
    for (std::shared_ptr<core::param::ParamSlot> param_slot : paramSlots) {
        Graph::ParamSlot psd;
        psd.class_name = std::string(param_slot->Name().PeekBuffer());
        psd.description = std::string(param_slot->Description().PeekBuffer());

        if (auto* p_ptr = param_slot->Param<core::param::ButtonParam>()) { psd.type = Graph::ParamType::BUTTON; }
        else if (auto* p_ptr = param_slot->Param<core::param::BoolParam>()) { psd.type = Graph::ParamType::BOOL; }
        else if (auto* p_ptr = param_slot->Param<core::param::ColorParam>()) { psd.type = Graph::ParamType::COLOR; }
        else if (auto* p_ptr = param_slot->Param<core::param::EnumParam>()) { psd.type = Graph::ParamType::ENUM; }
        else if (auto* p_ptr = param_slot->Param<core::param::FilePathParam>()) { psd.type = Graph::ParamType::FILEPATH; }
        else if (auto* p_ptr = param_slot->Param<core::param::FlexEnumParam>()) { psd.type = Graph::ParamType::FLEXENUM; }
        else if (auto* p_ptr = param_slot->Param<core::param::FloatParam>()) { psd.type = Graph::ParamType::FLOAT; }
        else if (auto* p_ptr = param_slot->Param<core::param::IntParam>()) { psd.type = Graph::ParamType::INT; }
        else if (auto* p_ptr = param_slot->Param<core::param::StringParam>()) { psd.type = Graph::ParamType::STRING; }
        else if (auto* p_ptr = param_slot->Param<core::param::TernaryParam>()) { psd.type = Graph::ParamType::TERNARY; }
        else if (auto* p_ptr = param_slot->Param<core::param::TransferFunctionParam>()) { psd.type = Graph::ParamType::TRANSFERFUNCTION; }
        else if (auto* p_ptr = param_slot->Param<core::param::Vector2fParam>()) { psd.type = Graph::ParamType::VECTOR2F; }
        else if (auto* p_ptr = param_slot->Param<core::param::Vector3fParam>()) { psd.type = Graph::ParamType::VECTOR3F; }
        else if (auto* p_ptr = param_slot->Param<core::param::Vector4fParam>()) { psd.type = Graph::ParamType::VECTOR4F; }
        else { psd.type = Graph::ParamType::UNKNOWN; }

        mod.param_slots.emplace_back(psd);
    }

    // CallerSlots
    for (std::shared_ptr<core::CallerSlot> caller_slot : callerSlots) {
        Graph::CallSlot csd;
        csd.name = std::string(caller_slot->Name().PeekBuffer());
        csd.description = std::string(caller_slot->Description().PeekBuffer());
        csd.compatible_call_idxs.clear();
        csd.connected_calls.clear();
        csd.parent_module = std::make_shared<Graph::Module>(mod);

        SIZE_T callCount = caller_slot->GetCompCallCount();
        for (SIZE_T i = 0; i < callCount; ++i) {
            std::string comp_call_class_name = std::string(caller_slot->GetCompCallClassName(i));
            auto calls_cnt = this->calls_list.size();
            for (size_t idx = 0; idx < calls_cnt; ++idx) {
                if (this->calls_list[idx].class_name == comp_call_class_name) {
                    csd.compatible_call_idxs.emplace_back(idx);
                }
            }
        }

        mod.caller_slots.emplace_back(csd);
    }

    // CalleeSlots
    for (std::shared_ptr<core::CalleeSlot> callee_slot : calleeSlots) {
        Graph::CallSlot csd;
        csd.name = std::string(callee_slot->Name().PeekBuffer());
        csd.description = std::string(callee_slot->Description().PeekBuffer());
        csd.compatible_call_idxs.clear();
        csd.connected_calls.clear();
        csd.parent_module = std::make_shared<Graph::Module>(mod);

        SIZE_T callbackCount = callee_slot->GetCallbackCount();
        std::vector<std::string> callNames, funcNames;
        std::set<std::string> uniqueCallNames, completeCallNames;
        for (SIZE_T i = 0; i < callbackCount; ++i) {
            uniqueCallNames.insert(callee_slot->GetCallbackCallName(i));
            callNames.push_back(callee_slot->GetCallbackCallName(i));
            funcNames.push_back(callee_slot->GetCallbackFuncName(i));
        }
        size_t ll = callNames.size();
        assert(ll == funcNames.size());
        for (std::string callName : uniqueCallNames) {
            core::factories::CallDescriptionManager::description_ptr_type desc = core_instance->GetCallDescriptionManager().Find(callName.c_str());
            bool allFound = true;
            if (desc != nullptr) {
                for (unsigned int i = 0; i < desc->FunctionCount(); ++i) {
                    std::string fn = desc->FunctionName(i);
                    bool found = false;
                    for (size_t j = 0; j < ll; ++j) {
                        if ((callNames[j] == callName) && (funcNames[j] == fn)) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        allFound = false;
                        break;
                    }
                }
            }
            else {
                allFound = false;
            }
            if (allFound) {
                completeCallNames.insert(callName);
            }
        }

        for (std::string callName : completeCallNames) {
            auto calls_cnt = this->calls_list.size();
            for (size_t idx = 0; idx < calls_cnt; ++idx) {
                if (this->calls_list[idx].class_name == callName) {
                    csd.compatible_call_idxs.emplace_back(idx);
                }
            }
        }

        mod.callee_slots.emplace_back(csd);
    }

    paramSlots.clear();
    callerSlots.clear();
    calleeSlots.clear();
    root_mod_ns->RemoveChild(new_mod);
    new_mod->SetAllCleanupMarks();
    new_mod->PerformCleanup();
    new_mod.reset();

    return true;
}


bool megamol::gui::Graph::read_call_data(Graph::Call& call_data, const std::shared_ptr<const megamol::core::factories::CallDescription> call_desc) {

    call_data.class_name = std::string(call_desc->ClassName());
    call_data.description = std::string(call_desc->Description());
    /// call_data.plugin_name is not available here (set above).
    call_data.connected_callee_slot.reset();
    call_data.connected_caller_slot.reset();

    return true;
}


bool megamol::gui::Graph::PROTOTYPE_SaveGraph(std::string project_filename, megamol::core::CoreInstance* core_instance) {

    if (core_instance == nullptr) {
        vislib::sys::Log::DefaultLog.WriteError("Pointer to CoreInstance is nullptr. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
        return false;
    }

    std::string confstr;
    // ------------------------------------------------------------------------
    {
        std::stringstream confInstances, confModules, confCalls, confParams;
        bool already_found_main_view = false;

        // Search for top most view
        for (auto& mod : this->modules_graph) {
            // Check for not connected calle_slots
            bool is_main_view = false;
            if (mod.is_view) {
                bool callee_connected = false;
                for (auto& cs : mod.callee_slots) {
                    for (auto& call : cs.connected_calls) {
                        if (call != nullptr) {
                            callee_connected = true;
                        }
                    }
                }
                if (!callee_connected) {
                    if (!already_found_main_view) {
                        confInstances << "mmCreateView(\"" << mod.instance << "\",\"" << mod.class_name << "\",\"" << mod.full_name << "\")\n";
                        already_found_main_view = true;
                        is_main_view = true;
                    }
                    else {
                        vislib::sys::Log::DefaultLog.WriteError("Only one main view is allowed. Found another one. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
                    }
                }
            }
            if (!is_main_view) {
                confModules << "mmCreateModule(\"" << mod.class_name << "\",\"" << mod.full_name << "\")\n";
            }

            /*
            for (auto& p : m.param_slots) {
                if (p.type != Graph::ParamType::BUTTON) {
                    std::string val = slot->Parameter()->ValueString().PeekBuffer();
                    // Encode to UTF-8 string
                    vislib::StringA valueString;
                    vislib::UTF8Encoder::Encode(valueString, vislib::StringA(val.c_str()));
                    val = valueString.PeekBuffer();
                    confParams << "mmSetParamValue(\"" << p.full_name << "\",[=[" << val << "]=])\n";
                }
            }
            */

            for (auto& cr : mod.caller_slots) {
                for (auto& call : cr.connected_calls) {
                    confCalls << "mmCreateCall(\"" << call->class_name << "\",\"" 
                        << mod.full_name  << "::" << cr.name << "\",\""
                        << call->connected_callee_slot->parent_module->full_name << "::" << call->connected_callee_slot->name << "\")\n";
                }
            }
        }
       
        confstr = confInstances.str() + "\n" + confModules.str() + "\n" + confCalls.str() + "\n" + confParams.str() + "\n";
    }
    // ------------------------------------------------------------------------

    try {
        std::ofstream file;
        file.open(project_filename);
        if (file.good()) {
            file << confstr.c_str();
            file.close();
        }
        else {
            vislib::sys::Log::DefaultLog.WriteError("Unable to create project file. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);
            file.close();
            return false;
        }
    }
    catch (...) {
    }

    return true;
}
