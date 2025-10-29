# Check if dsp-toolkit is available via pkg-config
ifneq ($(shell pkg-config --exists dsp-toolkit && echo yes),yes)
    $(error dsp-toolkit not found. Run 'nix develop' first)
endif

TOOLKIT_CFLAGS = $(shell pkg-config --cflags dsp-toolkit)
TOOLKIT_LIBDIR = $(shell pkg-config --variable=libdir dsp-toolkit)
TOOLKIT_LIBPATH = $(TOOLKIT_LIBDIR)/lib-DSP-Toolkit.so

# Helper function to run ROOT macros with toolkit loaded
define run_root_macro
	@cd macros && \
	 echo "Loading DSP Toolkit from: $(TOOLKIT_LIBDIR)" && \
	 echo "Executing $(1) macro..." && \
	 root -l -b -q -e 'gSystem->AddIncludePath("$(TOOLKIT_CFLAGS)"); \
	                    gSystem->Load("$(TOOLKIT_LIBPATH)");' \
			    $(2).C
endef

initial:
	$(call run_root_macro,initial processing,InitialProcessing)

calibrate:
	$(call run_root_macro,calibration,Calibration)

background:
	$(call run_root_macro,background subtraction,Background)

plots:
	$(call run_root_macro,main plotting,Plotting)

PSD:
	$(call run_root_macro,pulse shape discrimination,PSD)

optimize:
	$(call run_root_macro,PSD gate optimization,OptimizeGates)

clean:
	@echo "Nothing to clean (macros run interpreted)"

.PHONY: initial calibrate background plots PSD optimize clean

