#ifndef INTERACTIVEFITEDITOR_H
#define INTERACTIVEFITEDITOR_H

#include "FittingUtils.hpp"
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
#include <TLatex.h>
#include <TMath.h>
#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TSystem.h>
#include <TTimer.h>
#include <iostream>

/**
 * @brief Interactive editor for a `TF1` fit from FittingUtils.
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
class InteractiveFitEditor : public TGMainFrame {
private:
  static const Int_t kSliderRes = 10000;
  static const Int_t kNDrawPts = 500;

  static const Int_t kBtnRefit = 1000;
  static const Int_t kBtnAccept = 1001;
  static const Int_t kBtnCancel = 1002;
  static const Int_t kBtnReset = 1003;
  static const Int_t kSliderBase = 2000;
  static const Int_t kEntryBase = 3000;
  static const Int_t kFixBase = 4000;
  static const Int_t kLoBoundBase = 5000;
  static const Int_t kHiBoundBase = 6000;
  static const Int_t kRangeSlider = 7000;
  static const Int_t kRangeLoEntry = 7001;
  static const Int_t kRangeHiEntry = 7002;

  TH1 *hist_;
  TF1 *fit_func_;
  TString info_label_text_;
  Double_t range_low_;
  Double_t range_high_;
  Double_t original_range_low_;
  Double_t original_range_high_;
  Double_t hist_x_min_;
  Double_t hist_x_max_;
  Int_t num_peaks_;
  Int_t num_params_;

  Double_t *original_params_;
  Double_t *original_bounds_low_;
  Double_t *original_bounds_high_;
  Bool_t *original_fixed_;

  Double_t *current_bounds_low_;
  Double_t *current_bounds_high_;

  TRootEmbeddedCanvas *embedded_canvas_;
  TPad *main_pad_;
  TPad *residual_pad_;

  TGHSlider **sliders_;
  TGNumberEntry **value_entries_;
  TGCheckButton **fix_checks_;
  TGNumberEntry **lo_bound_entries_;
  TGNumberEntry **hi_bound_entries_;

  TGDoubleHSlider *range_slider_;
  TGNumberEntry *range_lo_entry_;
  TGNumberEntry *range_hi_entry_;

  TH1 *hist_draw_;
  TGraph *comp_graphs_[3][4];
  TF1 *bkg_draw_;
  TGraph *res_graph_;
  TF1 *zero_line_;
  TF1 *plus3_line_;
  TF1 *minus3_line_;
  TLatex *chi2_label_;
  Int_t n_res_points_;

  Bool_t needs_redraw_;
  Bool_t accepted_;
  Bool_t done_;
  Bool_t syncing_;
  TTimer *redraw_timer_;

  void BuildGUI();
  void BuildPeakTab(TGCompositeFrame *parent, Int_t peak_idx);
  void BuildBackgroundTab(TGCompositeFrame *parent);
  void AddParamRow(TGCompositeFrame *parent, Int_t param_idx, const char *name);

  void InitDrawing();
  void UpdateCanvas();
  void UpdateCompPoints();
  void UpdateResPoints();

  void SyncAllWidgets();
  void SyncWidget(Int_t param_idx);

  void OnSliderMoved(Int_t param_idx);
  void OnEntryChanged(Int_t param_idx);
  void OnBoundsChanged(Int_t param_idx);
  void OnFixToggled(Int_t param_idx);
  void OnRangeChanged();

  void DoRefit();
  void DoAccept();
  void DoCancel();
  void DoReset();

  Int_t ValToSlider(Int_t param_idx, Double_t val);
  Double_t SliderToVal(Int_t param_idx, Int_t pos);
  Bool_t IsFixed(Int_t param_idx);
  void GetDefaultBounds(Int_t param_idx, Double_t &lo, Double_t &hi);
  Int_t BkgConstIdx() { return num_peaks_ * 10; }
  Int_t BkgSlopeIdx() { return num_peaks_ * 10 + 1; }
  static Int_t PeakStyle(Int_t peak_idx);

public:
  /**
   * @brief Build the editor around an existing histogram and fit function.
   * @param parent     Parent window, normally `gClient->GetRoot()`.
   * @param hist       Histogram being fitted. Borrowed.
   * @param fit_func   Fit function, edited in place. Borrowed.
   * @param range_low  Initial lower fit bound.
   * @param range_high Initial upper fit bound.
   * @param num_peaks  Peaks in the model, 1 to 3; decides the tab layout.
   * @param info_label Optional annotation shown in the editor.
   */
  InteractiveFitEditor(const TGWindow *parent, TH1 *hist, TF1 *fit_func,
                       Double_t range_low, Double_t range_high, Int_t num_peaks,
                       const TString &info_label = "");
  /// @brief Destroys the widgets and drawing objects the editor created.
  virtual ~InteractiveFitEditor();

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
  /// @brief The canvas holding the fit and residual pads.
  /// @return Borrowed pointer; the editor owns it.
  TRootEmbeddedCanvas *GetEmbeddedCanvas() { return embedded_canvas_; }
};

/**
 * @brief Open the editor and pump its event loop until the user is done.
 *
 * Disables ROOT batch mode for the duration and restores it afterwards, so
 * ordinary plot output is unaffected, and installs the tolerant X error handler
 * from InteractiveEditorX11Guard.hpp around the loop.
 *
 * @param hist       Histogram being fitted.
 * @param fit_func   Fit function, updated in place if the user accepts.
 * @param range_low  Initial lower fit bound.
 * @param range_high Initial upper fit bound.
 * @param num_peaks  Peaks in the model, 1 to 3.
 * @param info_label Optional annotation shown in the editor.
 *
 * @return `kTRUE` if the user accepted; `kFALSE` on cancel, in which case
 *         @p fit_func is left as it was.
 */
Bool_t LaunchInteractiveFitEditor(TH1 *hist, TF1 *fit_func, Double_t range_low,
                                  Double_t range_high, Int_t num_peaks,
                                  const TString &info_label);

#endif
