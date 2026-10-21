# Void run — 2026-09-22

This directory holds the first attempt at the rog prospective test. **It does
not satisfy the registered protocol and must not be analysed.**

    registered:  20 backbones x 4 ProteinMPNN draws = 80 folds
    actually ran: 16 backbones, 1-2 draws each      = 20 folds

Cause: every batch's ProteinMPNN output holds one sequence per backbone. The
prepare script took whatever existed instead of failing. The same trap was
already worked around in `root_breadth.sh` and the workaround was not carried
over here.

Two consequences beyond the fold count. Best-of-1 is far noisier than
best-of-4, which attenuates the correlation under test — so a real effect
would be pushed toward null. And four backbones dropped out entirely, so the
registered group composition was not what ran.

The interim read from this run (compact 0.169 vs extended 0.130, gap +0.039,
p = 0.064) was reported and then retracted. It is **not** evidence that rog
failed to replicate; the test had not been performed.

`REGISTRATION.md` in the parent directory is unedited and still governs.
