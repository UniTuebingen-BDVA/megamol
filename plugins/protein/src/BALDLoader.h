/**
 * @author Nils Lichtenberg
 */
#ifndef MMPROTEINPLUGIN_BALDLOADER_H_INCLUDED
#define MMPROTEINPLUGIN_BALDLOADER_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/param/ParamSlot.h"
#include "mmstd/data/AnimDataModule.h"
#include <fstream>
#include <utility>

namespace megamol::protein {

class BALDLoader : public megamol::core::view::AnimDataModule {
public:
    /** Ctor */
    BALDLoader();

    /** Dtor */
    ~BALDLoader() override;

    /**
     * Answer the name of this module.
     *
     * @return The name of this module.
     */
    static const char* ClassName() {
        return "BALDLoader";
    }

    /**
     * Answer a human readable description of this module.
     *
     * @return A human readable description of this module.
     */
    static const char* Description() {
        return "Loads binary atom layer data (BALD).";
    }

    /**
     * Answers whether this module is available on the current system.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() {
        return true;
    }

protected:
    /**
     * Implementation of 'Create'.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    bool create() override;

    /**
     * Call callback to get the data
     *
     * @param c The calling call
     *
     * @return True on success
     */
    bool getData(core::Call& call);

    /**
     * Call callback to get the extent of the data
     *
     * @param c The calling call
     *
     * @return True on success
     */
    bool getExtent(core::Call& call);

    /**
     * Call callback to check whether data has been changed/needs update
     *
     * @param c The calling call
     *
     * @return whether data gas changed
     */
    bool dataChanged(core::Call& call) {
        return false; /*return solventResidues.IsDirty();*/
    }

    /**
     * Implementation of 'Release'.
     */
    void release() override;

    /**
     * Creates a frame to be used in the frame cache. This method will be
     * called from within 'initFrameCache'.
     *
     * @return The newly created frame object.
     */
    Frame* constructFrame() const override;

    /**
     * Loads one frame of the data set into the given 'frame' object. This
     * method may be invoked from another thread. You must take
     * precausions in case you need synchronised access to shared
     * ressources.
     *
     * @param frame The frame to be loaded.
     * @param idx The index of the frame to be loaded.
     */
    void loadFrame(Frame* frame, unsigned int idx) override;

    /** The pdb file name slot */
    core::param::ParamSlot baldFilenameSlot;

    /** The slot for requesting data */
    core::CalleeSlot getDataSlot;

    /**
     * Data provided by a .BALD file header.
     */
    struct HeaderData {
        int magicNumber, numAtoms, numFrames;
        float sasRadius;
    } headerData{};

    /**
     * Storage of frame data
     */
    class Frame : public megamol::core::view::AnimDataModule::Frame {
    public:
        /** Ctor */
        Frame(megamol::core::view::AnimDataModule& owner);

        /** Dtor */
        ~Frame() override;

        /**
         *
         * Set the frame Index.
         *
         * @param idx the index
         */
        void setFrameIdx(unsigned int idx) {
            this->frame = idx;
        }

        std::vector<float> getLayers() const {
            return layers;
        }
        void setLayers(std::vector<float> layersIn) {
            Frame::layers = std::move(layersIn);
        }
        unsigned int getMaxLayer() const {
            return maxLayer;
        }
        void setMaxLayer(unsigned int maxLayerIn) {
            Frame::maxLayer = maxLayerIn;
        }

    private:
        std::vector<float> layers;
        unsigned int maxLayer{};
    };
    bool readHeaderData();
};


} // namespace megamol::protein

#endif // MMPROTEINPLUGIN_BALDLOADER_H_INCLUDED
