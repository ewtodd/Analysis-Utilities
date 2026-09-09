BinaryUtils {#binary-utils}
===========

Read binary data files from CAEN acquisition software and the SOLARIS DAQ.

[TOC]

## CoMPASSReader

Reads CoMPASS binary files (`.bin`):

- Parses event headers, timestamps, energy values, and waveforms
- Supports all CoMPASS waveform codes (INPUT, RC_CR, TRAPEZOID, CFD, etc.)
- Decodes status flags (pileup, saturation, deadtime, trigger lost, etc.)

## WaveDump742Reader

Reads WaveDump binary files for DT5742 family digitizers:

- Parses event structure including DC offset and start index cell
- Optional timing corrections support

## SOLReader

Reads SOL binary files produced by the
[SOLARIS DAQ](https://github.com/goluckyryan/SOLARIS_DAQ):

- Parses all SOL data types (ALL, OneTrace, NoTrace, Minimum, MiniWithFineTime,
  Raw)
- Per-block access to energy, fine timestamp, flags, and probe types, plus
  `getAnalog0` / `getAnalog1`, `getDigital`, and `getOneTrace` accessors for
  trace-carrying formats
- `SetSkipTraces()` skips the trace payload (blocks are still parsed); blocks
  with an implausible trace length (> 100k samples) are skipped as a
  memory-safety measure
- `ToHit()` reduces the current block to a lightweight `SOLHit` (header fields
  only, no trace)
- Static `SplitSolFileByTime(input, output_dir, chunk_seconds)` splits a run
  file into time-bounded chunks named `<basename>_chunk<NNN>.sol`; reads and
  writes are bulk-buffered (4 MiB) so splitting runs close to disk speed
