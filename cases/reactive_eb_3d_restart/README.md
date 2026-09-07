# Reactive planar EB 3D checkpoint/restart

This case qualifies a serial, single-level, axis-aligned planar embedded
boundary continuation with chemistry and molecular transport enabled.  Run
`reference.nml` continuously, run `checkpoint_stop.nml` to stop after the
first committed step, then run `restart.nml` in the same directory.

The regression gate requires the restarted final CSV to be byte-identical to
the uninterrupted reference.  The stopped CSV must represent an earlier time,
and the checkpoint must contain its terminal marker.
