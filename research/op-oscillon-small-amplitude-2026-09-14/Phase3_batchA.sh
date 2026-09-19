#!/bin/sh
# Etap A (FROZEN): 13 biegow h=0.05 do t=2000, rownolegle (32 CPU).
# Uruchamiane DOPIERO po PASS Phase 2 (LOCK sec. 4).
B="C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/research/op-oscillon-small-amplitude-2026-09-14"
for a in 0.02 0.05 0.08 0.10; do
  for s in 3 6 10; do
    python "$B/Phase3_evolve.py" run "g_a${a}_s${s}" h05 2000 >> "$B/Phase3_batchA.log" 2>&1 &
  done
done
python "$B/Phase3_evolve.py" run vac h05 2000 >> "$B/Phase3_batchA.log" 2>&1 &
wait
echo "BATCH A COMPLETE"
