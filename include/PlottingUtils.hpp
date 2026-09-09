#ifndef PLOTTINGUTILS_H
#define PLOTTINGUTILS_H

#include <TCanvas.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

/// @brief Which axis scalings SaveFigure() should write out.
enum class PlotSaveOptions {
  kLINEAR, ///< Linear y only.
  kLOG,    ///< Logarithmic y only, written with a `log_` filename prefix.
  kBOTH    ///< Both, as two separate files.
};

/// @brief Output file format for saved figures.
enum class PlotSaveFormat {
  kPNG, ///< PNG raster output; line width 2.
  kPDF  ///< PDF vector output; line width 1.
};

/**
 * @brief All-static helpers for publication-quality ROOT graphics.
 *
 * No instantiation is required or possible in practice — every member is
 * static. SetStylePreferences() must be called before anything else; the other
 * methods print a warning to stdout if it has not been, but still run.
 *
 * @note Several methods return heap-allocated ROOT objects. ROOT's own
 *       ownership rules apply and are called out per method; in general the
 *       objects returned here are owned by the caller or by the current pad,
 *       never by this class.
 */
class PlottingUtils {
public:
  /**
   * @brief Install the global ROOT style and choose the output format.
   *
   * Sets `gStyle` globally: no stat or fit box, 0.06 axis title and label
   * sizes, 1.2 title offsets, grids on in a light grey dashed style, ticks on
   * both axes, and a 255-colour turbo palette for 2-D colour scales. Also fixes
   * the line width used by the Configure* methods — 2 for PNG, 1 for PDF.
   *
   * @param save_format Output format for SaveFigure(). Defaults to
   *                    PlotSaveFormat::kPNG.
   *
   * @note Mutates ROOT global state (`gStyle`) for the whole process.
   *       InitUtils::SetROOTPreferences() calls this and is the recommended
   *       entry point.
   */
  static void
  SetStylePreferences(PlotSaveFormat save_format = PlotSaveFormat::kPNG);

  /// @brief ConfigureGraph() followed by `Draw()`. No-op if @p graph is null.
  /// @param graph Graph to configure and draw onto the current pad.
  /// @param color ROOT colour index for the line.
  /// @param title Graph title, in the usual ROOT `"title;x;y"` form.
  static void ConfigureAndDrawGraph(TGraph *graph, Int_t color,
                                    const TString title = "");

  /// @brief ConfigureHistogram() followed by `Draw("HIST")`. No-op if null.
  /// @param hist  Histogram to configure and draw onto the current pad.
  /// @param color ROOT colour index for the line.
  /// @param title Histogram title, in the usual ROOT `"title;x;y"` form.
  static void ConfigureAndDrawHistogram(TH1 *hist, Int_t color,
                                        const TString title = "");

  /// @brief Configure2DHistogram() followed by `Draw("COLZ")`. No-op if null.
  /// @param hist   2-D histogram to configure and draw.
  /// @param canvas Canvas whose log-z and right margin are adjusted.
  /// @param title  Histogram title, in the usual ROOT `"title;x;y;z"` form.
  static void ConfigureAndDraw2DHistogram(TH2 *hist, TCanvas *canvas,
                                          const TString title = "");

  /**
   * @brief Apply the house style to a graph without drawing it.
   *
   * Sets line colour and width, 0.06 axis title and label sizes, 1.2 title
   * offsets, and 506 x-axis divisions.
   *
   * @param graph Graph to configure. Must not be null.
   * @param color ROOT colour index for the line.
   * @param title Graph title, in the usual ROOT `"title;x;y"` form.
   */
  static void ConfigureGraph(TGraph *graph, Int_t color,
                             const TString title = "");

  /**
   * @brief Apply the house style to a graph with error bars.
   *
   * As the `TGraph` overload, plus marker style 20 at size 1.2 in the same
   * colour as the line — the usual look for a data series with uncertainties.
   *
   * @param graph Graph to configure. Must not be null.
   * @param color ROOT colour index for both line and markers.
   * @param title Graph title, in the usual ROOT `"title;x;y"` form.
   */
  static void ConfigureGraph(TGraphErrors *graph, Int_t color,
                             const TString title = "");

  /**
   * @brief Apply the house style to a 1-D histogram without drawing it.
   *
   * Sets line colour and width, hollow fill, 0.06 axis title and label sizes,
   * 1.2 title offsets, and axis divisions tuned for spectra. The x axis is
   * forced out of exponent notation.
   *
   * @param hist  Histogram to configure. Null is tolerated and ignored.
   * @param color ROOT colour index for the line.
   * @param title Histogram title, in the usual ROOT `"title;x;y"` form.
   */
  static void ConfigureHistogram(TH1 *hist, Int_t color,
                                 const TString title = "");

  /**
   * @brief Apply the house style to a 2-D histogram and its canvas.
   *
   * Styles the axes as for the 1-D case, then enables log-z on @p canvas and
   * widens its right margin to 0.15 to make room for the colour axis.
   *
   * @param hist   Histogram to configure. Null is tolerated and ignored.
   * @param canvas Canvas to adjust. Null is tolerated and ignored.
   * @param title  Histogram title, in the usual ROOT `"title;x;y;z"` form.
   */
  static void Configure2DHistogram(TH2 *hist, TCanvas *canvas,
                                   const TString title = "");

  /**
   * @brief Create a 1200x800 canvas with grid and ticks already set up.
   *
   * @param logy `kTRUE` to start with a logarithmic y axis.
   *
   * @return A newly allocated `TCanvas` with a randomised name (see
   *         GetRandomName()). **The caller owns this pointer.**
   *
   * @note The canvas is made current, so subsequent `Draw()` calls land on it.
   */
  static TCanvas *GetConfiguredCanvas(Bool_t logy = kFALSE);

  /**
   * @brief Write a canvas to disk under the configured plots base directory.
   *
   * The file lands in `<plots_base>/` or `<plots_base>/<subdirectory>/`, with
   * the extension chosen by SetStylePreferences(). Missing directories are
   * created. For PlotSaveOptions::kLOG and PlotSaveOptions::kBOTH the log-y
   * variant is written with a `log_` filename prefix.
   *
   * @param canvas              Canvas to save. Must not be null.
   * @param output_filename     Base filename, without extension.
   * @param output_subdirectory Subdirectory under the plots base; empty means
   *                            the base itself.
   * @param save_options        Which axis scalings to write.
   *                            PlotSaveOptions::kBOTH by default.
   *
   * @warning Requesting a logarithmic variant for a canvas containing a `TH2`
   *          is treated as a usage error and **terminates the process** via
   *          `std::exit(1)` after printing an explanation. The log option
   *          refers to the y axis and is meaningless for a 2-D colour plot,
   *          whose z axis is already logarithmic if it went through
   *          Configure2DHistogram(). Use PlotSaveOptions::kLINEAR for those.
   *
   * @note The canvas is left with log-y disabled on return, whatever it was
   *       set to beforehand.
   */
  static void SaveFigure(TCanvas *canvas, TString output_filename,
                         TString output_subdirectory = "",
                         PlotSaveOptions save_options = PlotSaveOptions::kBOTH);

  /**
   * @brief Set the base directory that saved figures are written under.
   *
   * @param dir Base directory. Pass an absolute path so output is anchored to a
   *            project root regardless of the current working directory.
   *            Trailing slashes are stripped.
   *
   * @note The default is the CWD-relative `"plots"`.
   *       InitUtils::SetROOTPreferences() calls this for you.
   */
  static void SetPlotsBaseDir(const TString &dir);

  /// @brief Current base directory for saved figures, without trailing slash.
  /// @return The directory last set by SetPlotsBaseDir(), or `"plots"`.
  static TString GetPlotsBaseDir();

  /**
   * @brief Draw a fit over its data with a residual panel and save it.
   *
   * Produces a two-pad figure: the histogram with the total fit and each
   * component overlaid on top (70% of the height), and a pull panel beneath
   * (30%) marked with dashed guide lines at ±3σ. The displayed x range is
   * padded 10% either side of the fit range.
   *
   * @param hist                Data histogram. Must not be null.
   * @param total_graph         Curve for the summed model.
   * @param component_graphs    One curve per model component; may be empty.
   * @param fit_range_low       Lower fit bound, in the histogram's x units.
   * @param fit_range_high      Upper fit bound, in the histogram's x units.
   * @param output_name         Base filename, without extension.
   * @param output_subdirectory Subdirectory under the plots base. `"fits"` by
   *                            default.
   * @param label               Optional annotation drawn on the plot; empty
   *                            draws nothing.
   * @param logy                `kTRUE` (default) for a logarithmic y axis on
   *                            the upper pad.
   *
   * @note Saves through SaveFigure(), so the same base-directory and format
   *       rules apply. A companion pull histogram is written alongside.
   */
  static void
  PlotFitWithResiduals(TH1 *hist, TGraph *total_graph,
                       const std::vector<TGraph *> &component_graphs,
                       Float_t fit_range_low, Float_t fit_range_high,
                       const TString &output_name,
                       const TString &output_subdirectory = "fits",
                       const TString &label = "", Bool_t logy = kTRUE);

  /**
   * @brief Create, style and draw a legend on the current pad.
   *
   * @param x1 Left edge, in NDC (0-1).
   * @param x2 Right edge, in NDC (0-1).
   * @param y1 Bottom edge, in NDC (0-1).
   * @param y2 Top edge, in NDC (0-1).
   *
   * @return The drawn `TLegend`, so entries can be added to it. Owned by the
   *         current pad once drawn.
   *
   * @warning The argument order is `(x1, x2, y1, y2)` — both x bounds, then
   *          both y bounds — **not** ROOT's own `TLegend(x1, y1, x2, y2)`.
   */
  static TLegend *AddLegend(Double_t x1 = 0.7, Double_t x2 = 0.9,
                            Double_t y1 = 0.7, Double_t y2 = 0.9);

  /**
   * @brief Create and draw a text annotation on the current pad.
   *
   * Useful for subplot labels such as "(a)" and "(b)", or any free annotation.
   * The text is right-top aligned at the given point.
   *
   * @param label Text to draw. ROOT `TLatex` markup is honoured.
   * @param x     Horizontal position, in NDC (0-1).
   * @param y     Vertical position, in NDC (0-1).
   * @param angle Rotation in degrees, counter-clockwise. 0 by default.
   *
   * @return The drawn `TLatex`. Owned by the current pad once drawn.
   */
  static TLatex *AddText(const TString label, Double_t x = 0.9,
                         Double_t y = 0.85, Double_t angle = 0);

  /// @brief A palette of visually distinct ROOT colour indices.
  /// @return 24 colour indices, intended to be indexed by series number.
  static std::vector<Int_t> GetDefaultColors();

  /**
   * @brief Generate a name unlikely to collide with existing ROOT objects.
   *
   * ROOT keeps named objects in global lists and warns (or misbehaves) on
   * duplicate names; this sidesteps that for throwaway canvases.
   *
   * @return A name of the form `nameX.XXXXXXX` from a clock-seeded generator.
   */
  static TString GetRandomName();

  /// @brief Line width the Configure* methods apply.
  /// @return 2 after SetStylePreferences(PlotSaveFormat::kPNG), 1 for kPDF.
  static Width_t GetLineWidth() { return line_width_; };

private:
  static PlotSaveFormat save_format_;
  static Bool_t preferences_set_;
  static Width_t line_width_;
  static TString plots_base_dir_;
  static void SetTurboPalette();
  static void WarnIfNotConfigured(const TString method_name);
  static void PlotPullHistogram(TGraph *residuals, const TString &output_name,
                                const TString &output_subdirectory);
};

#endif
