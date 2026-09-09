#ifndef INTERACTIVESIMULTANEOUSFITEDITOR_H
#define INTERACTIVESIMULTANEOUSFITEDITOR_H

#include "RooFitUtils.hpp"

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <TCanvas.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TGDoubleSlider.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGSlider.h>
#include <TGTab.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLatex.h>
#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TTimer.h>
#include <map>
#include <vector>

/**
 * @brief Everything the simultaneous editor needs to show one channel.
 *
 * The caller fills in the model and data members; the editor fills in the
 * widget and drawing members as it builds that channel's tab, and writes edited
 * values back through @ref params on accept.
 *
 * @note Every pointer here is borrowed. RooFitUtils owns the model objects and
 *       the editor owns the widgets, so a view must not outlive either.
 */
struct SimEditorChannelView {
  TString name; ///< Channel name, shown on the tab.
  TH1 *hist;    ///< Display histogram for this channel.
  /// Event-level values behind @ref hist, so the display can be rebinned live
  /// as the fit range changes.
  const std::vector<Double_t> *events;
  Float_t display_bin_width_kev;       ///< Display bin width.
  RooAbsPdf *pdf;                      ///< This channel's summed model.
  RooAbsData *data;                    ///< This channel's unbinned dataset.
  std::vector<RooFitPeakModel> *peaks; ///< Per-peak parameter models.
  RooFitBackgroundModel *bkg;          ///< Background model.
  Int_t num_peaks;                     ///< Peaks in this channel.

  /// @name Parameter state
  /// @ref params is the live set the editor writes through; the `original_*`
  /// vectors are the snapshot Reset restores, and the `current_*` vectors track
  /// bounds as the user edits them.
  /// @{
  std::vector<RooRealVar *> params;
  std::vector<Double_t> original_params;
  std::vector<Double_t> original_bounds_low;
  std::vector<Double_t> original_bounds_high;
  std::vector<Bool_t> original_fixed;
  std::vector<Double_t> current_bounds_low;
  std::vector<Double_t> current_bounds_high;
  /// @}

  /// @name Widgets and drawing objects
  /// Populated by the editor when it builds this channel's tab.
  /// @{

  std::vector<TGHSlider *> sliders;
  std::vector<TGNumberEntry *> value_entries;
  std::vector<TGCheckButton *> fix_checks;
  std::vector<TGNumberEntry *> lo_bound_entries;
  std::vector<TGNumberEntry *> hi_bound_entries;

  TRootEmbeddedCanvas *embedded_canvas;
  TPad *main_pad;
  TPad *residual_pad;
  TH1 *hist_draw;
  TGraph *total_graph;
  TGraph *bkg_graph;
  TGraph *comp_graphs[3][4];
  TGraph *res_graph;
  TF1 *zero_line;
  TF1 *plus3_line;
  TF1 *minus3_line;
  TLatex *chi2_label;
  Int_t n_res_points;
  /// @}
};

/**
 * @brief Interactive editor for a multi-channel simultaneous RooFit fit.
 *
 * The fit and residual pads update live as parameters move. Parameters can be
 * fixed or freed individually, bounds edited, and Refit re-runs the minimiser
 * from the current values as starting points. The residual panel is marked with
 * dashed guide lines at plus and minus three sigma.
 *
 * @note Normally reached through the Launch function below rather than
 *       constructed directly; the launcher owns the event loop, the batch-mode
 *       toggle and the X error-handler guard.
 *
 * @warning None of the objects passed in are owned by the editor, and all of
 *          them must outlive it. Accepting writes the edited values back into
 *          them in place.
 */
class InteractiveSimultaneousFitEditor : public TGMainFrame {
private:
  static const Int_t kSliderRes = 10000;
  static const Int_t kNDrawPts = 500;

  static const Int_t kBtnRefit = 1000;
  static const Int_t kBtnAccept = 1001;
  static const Int_t kBtnCancel = 1002;
  static const Int_t kBtnReset = 1003;

  static const Int_t kRangeSlider = 7000;
  static const Int_t kRangeLoEntry = 7001;
  static const Int_t kRangeHiEntry = 7002;

  static const Int_t kSliderStride = 100;
  static const Int_t kSliderBase = 10000;
  static const Int_t kEntryBase = 20000;
  static const Int_t kFixBase = 30000;
  static const Int_t kLoBoundBase = 40000;
  static const Int_t kHiBoundBase = 50000;

  RooSimultaneous *sim_pdf_;
  RooAbsData *combined_data_;
  RooRealVar *x_;
  std::vector<SimEditorChannelView> channels_;
  TString info_label_text_;
  Double_t range_low_;
  Double_t range_high_;
  Double_t original_range_low_;
  Double_t original_range_high_;
  Double_t hist_x_min_;
  Double_t hist_x_max_;

  TGDoubleHSlider *range_slider_;
  TGNumberEntry *range_lo_entry_;
  TGNumberEntry *range_hi_entry_;

  Bool_t needs_redraw_;
  Bool_t accepted_;
  Bool_t done_;
  Bool_t syncing_;
  Bool_t fit_debug_;
  TTimer *redraw_timer_;

  Int_t ChannelIndexFromWidgetId(Int_t base, Int_t parm1, Int_t &local_idx);

  void BuildGUI();
  void BuildChannelTab(TGCompositeFrame *parent, Int_t ch_idx);
  void BuildPeakSubTab(TGCompositeFrame *parent, Int_t ch_idx, Int_t peak_idx);
  void BuildBackgroundSubTab(TGCompositeFrame *parent, Int_t ch_idx);
  void AddParamRow(TGCompositeFrame *parent, Int_t ch_idx, Int_t param_idx,
                   const char *name);

  void InitDrawing();
  void InitChannelDrawing(SimEditorChannelView &cv);
  void UpdateCanvases();
  void UpdateChannelGraphs(SimEditorChannelView &cv);
  void UpdateChannelResiduals(SimEditorChannelView &cv);
  void UpdateChannelChi2(SimEditorChannelView &cv);

  void SyncAllWidgets();
  void SyncChannelWidget(Int_t ch_idx, Int_t param_idx);

  void OnSliderMoved(Int_t ch_idx, Int_t param_idx);
  void OnEntryChanged(Int_t ch_idx, Int_t param_idx);
  void OnBoundsChanged(Int_t ch_idx, Int_t param_idx);
  void OnFixToggled(Int_t ch_idx, Int_t param_idx);
  void OnRangeChanged();

  void DoRefit();
  void DoAccept();
  void DoCancel();
  void DoReset();

  // Re-tie each channel's linear-background slope lower bound to the current
  // shared fit range so 1+slope*x stays positive across the whole window.
  void ApplyBackgroundSlopeBounds();

  // Debug aid: scan each channel's components for a non-finite normalization
  // integral or a non-finite raw value across the range, to name the pdf
  // behind an invalid NLL. Only called when fit_debug_ is set.
  void DiagnoseInvalidComponents();

  Int_t ValToSlider(Int_t ch_idx, Int_t param_idx, Double_t val);
  Double_t SliderToVal(Int_t ch_idx, Int_t param_idx, Int_t pos);
  Bool_t IsFixed(Int_t ch_idx, Int_t param_idx);
  void GetDefaultBounds(Int_t ch_idx, Int_t param_idx, Double_t &lo,
                        Double_t &hi);
  Int_t BkgConstIdx(Int_t ch_idx) { return channels_[ch_idx].num_peaks * 10; }
  Int_t BkgSlopeIdx(Int_t ch_idx) {
    return channels_[ch_idx].num_peaks * 10 + 1;
  }
  static Int_t PeakStyle(Int_t peak_idx);

public:
  /**
   * @brief Build the editor around a converged simultaneous model.
   * @param parent         Parent window, normally `gClient->GetRoot()`.
   * @param sim_pdf        The `RooSimultaneous` model. Borrowed.
   * @param combined_data  Combined dataset across every channel. Borrowed.
   * @param x              Shared observable. Borrowed.
   * @param channel_views  One view per channel; each gets its own tab, and its
   *                       parameters are edited in place.
   * @param range_low      Initial lower fit bound.
   * @param range_high     Initial upper fit bound.
   * @param info_label     Optional annotation shown in the editor.
   * @param fit_debug      `kTRUE` to un-suppress RooFit evaluation errors and
   *                       print the seed NLL on refit, so an invalid-NLL
   *                       failure names the offending PDF.
   */
  InteractiveSimultaneousFitEditor(
      const TGWindow *parent, RooSimultaneous *sim_pdf,
      RooAbsData *combined_data, RooRealVar *x,
      const std::vector<SimEditorChannelView> &channel_views,
      Double_t range_low, Double_t range_high, const TString &info_label = "",
      Bool_t fit_debug = kFALSE);
  /// @brief Destroys the widgets and drawing objects the editor created.
  virtual ~InteractiveSimultaneousFitEditor();

  /**
   * @brief ROOT GUI message dispatch for every widget in the editor.
   * @param msg   Encoded message type and subtype.
   * @param parm1 Widget id that raised it.
   * @param parm2 Message-specific payload.
   * @return `kTRUE` once handled.
   */
  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
  /**
   * @brief Redraw tick.
   *
   * Parameter edits set a dirty flag rather than redrawing inline; this
   * coalesces them so dragging a slider does not queue one full redraw per
   * pixel of travel.
   *
   * @param timer Timer that fired.
   * @return `kTRUE` once handled.
   */
  virtual Bool_t HandleTimer(TTimer *timer);
  /// @brief Window-manager close. Treated as a cancel, not an accept.
  virtual void CloseWindow();

  /// @brief Whether the user accepted rather than cancelled.
  /// @return `kTRUE` if Accept was pressed. Only meaningful once IsDone().
  Bool_t WasAccepted() const { return accepted_; }
  /// @brief Whether the editor has finished and the loop may exit.
  Bool_t IsDone() const { return done_; }
  /// @brief The coalescing redraw timer, for the driving event loop.
  /// @return Borrowed pointer; the editor owns it.
  TTimer *GetRedrawTimer() { return redraw_timer_; }
  /// @brief The per-channel views, carrying the edited parameters.
  /// @return Borrowed pointer to the editor's own vector.
  std::vector<SimEditorChannelView> *GetChannels() { return &channels_; }
};

/**
 * @brief Open the simultaneous editor and pump its event loop.
 *
 * Same batch-mode and X error-handler handling as
 * LaunchInteractiveFitEditor(). Refit re-runs the joint minimisation across
 * every channel at once, not one channel at a time.
 *
 * @param sim_pdf       The `RooSimultaneous` model.
 * @param combined_data Combined dataset across every channel.
 * @param x             Shared observable.
 * @param channel_views One view per channel, updated on accept.
 * @param range_low     Initial lower fit bound.
 * @param range_high    Initial upper fit bound.
 * @param info_label    Optional annotation shown in the editor.
 * @param fit_debug     `kTRUE` to un-suppress RooFit evaluation errors.
 *
 * @return `kTRUE` if the user accepted; `kFALSE` on cancel.
 */
Bool_t LaunchInteractiveSimultaneousFitEditor(
    RooSimultaneous *sim_pdf, RooAbsData *combined_data, RooRealVar *x,
    std::vector<SimEditorChannelView> &channel_views, Double_t range_low,
    Double_t range_high, const TString &info_label, Bool_t fit_debug = kFALSE);

#endif
