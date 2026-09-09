#ifndef INTERACTIVEEDITORX11GUARD_H
#define INTERACTIVEEDITORX11GUARD_H

#include <Rtypes.h>
#include <dlfcn.h>

/**
 * @file InteractiveEditorX11Guard.hpp
 * @brief Stops a recoverable X protocol error from becoming a fatal SIGSEGV.
 *
 * When a `TGTextEntry` (the editors' number entries, for instance) is redrawn
 * through the TTF backend, `TGX11TTF::DrawString` issues an `XGetGeometry`
 * round-trip on the widget's drawable. If that drawable is momentarily stale
 * the X server answers with `BadWindow` / `BadDrawable`. Xlib then calls ROOT's
 * `RootX11ErrorHandler`, which looks the offending window up by id and calls a
 * virtual method on it to format the message — and if that window object is
 * already gone, the handler itself crashes. A recoverable protocol error thus
 * becomes a fatal SIGSEGV that takes down the entire ROOT session.
 *
 * While an interactive editor pumps its own event loop, these helpers install a
 * no-op handler that swallows such errors — Xlib ignores the return value for
 * non-fatal errors and carries on — and then restore the previous one.
 *
 * `XSetErrorHandler` is resolved at runtime via `dlsym`, so this library keeps
 * no link-time dependency on libX11; the symbol comes from libX11, which ROOT's
 * GUI libraries load transitively whenever a `gClient` exists.
 */

/// @brief Xlib error handler signature, with opaque pointers so that no libX11
///        headers are needed here.
typedef int (*AUXErrorHandler)(void *display, void *event);

/// @brief Signature of `XSetErrorHandler` itself, resolved via `dlsym`.
typedef AUXErrorHandler (*AUXSetErrorHandlerFn)(AUXErrorHandler);

/// @brief The tolerant handler: ignores the error and reports success.
/// @return Always `0`; Xlib disregards this for non-fatal errors.
inline int AUIgnoreXError(void *, void *) { return 0; }

/**
 * @brief What AUInstallTolerantXErrorHandler() needs in order to undo itself.
 *
 * Treat as opaque: pass it back to AURestoreXErrorHandler() unchanged.
 */
struct AUXErrorHandlerSave {
  /// Resolved `XSetErrorHandler`, or null if it could not be found.
  AUXSetErrorHandlerFn set_fn;
  /// Handler that was installed before, to be put back on restore.
  AUXErrorHandler prev;
};

/**
 * @brief Install the tolerant X error handler for the duration of an event
 * loop.
 *
 * Resolves `XSetErrorHandler` from the default symbol scope; if libX11 was
 * loaded with local scope the lookup is retried against an `RTLD_NOLOAD`
 * handle, which promotes the already-mapped library rather than loading a
 * second copy.
 *
 * @return State to hand to AURestoreXErrorHandler(). If the symbol could not be
 *         resolved the call is a silent no-op and restoring is harmless, so the
 *         result never needs checking.
 *
 * @warning This suppresses **all** non-fatal X protocol errors for the whole
 *          process while installed, not just the redraw race it targets. Keep
 *          the window as narrow as the editor loop, and always pair it with
 *          AURestoreXErrorHandler().
 */
inline AUXErrorHandlerSave AUInstallTolerantXErrorHandler() {
  AUXErrorHandlerSave s;
  s.set_fn = (AUXSetErrorHandlerFn)dlsym(RTLD_DEFAULT, "XSetErrorHandler");
  if (!s.set_fn) {
    // libX11 may have been brought in with local scope; promote and retry.
    void *h = dlopen("libX11.so.6", RTLD_NOW | RTLD_GLOBAL | RTLD_NOLOAD);
    if (h) {
      s.set_fn = (AUXSetErrorHandlerFn)dlsym(h, "XSetErrorHandler");
    }
  }
  s.prev = s.set_fn ? s.set_fn(AUIgnoreXError) : (AUXErrorHandler)0;
  return s;
}

/**
 * @brief Put back the handler that was in place before the tolerant one.
 *
 * @param s Value returned by AUInstallTolerantXErrorHandler(). Safe to pass a
 *          state whose resolution failed; the call then does nothing.
 */
inline void AURestoreXErrorHandler(AUXErrorHandlerSave s) {
  if (s.set_fn) {
    s.set_fn(s.prev);
  }
}

#endif
