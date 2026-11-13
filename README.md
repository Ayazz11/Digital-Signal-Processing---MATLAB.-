# Discrete-Time Signal Processing from Scratch (MATLAB)

Design-level implementations of core DSP building blocks, written from the mathematical definitions **without** calling `conv`, `cconv`, `fft`, `ifft` or `filter`. Built-in functions are used only as a *reference* to check the results.

## Status

| Module | Topic |
|---|---|---|
| 1 | Sequence generation, linear and circular convolution |
| 2 | N-point DFT / IDFT, zero-padding, truncation vs folding, DFT–DTFT relationship |
| 3 | Radix-2 FFT, fast convolution |
| 4 | Windowed-sinc FIR design, Butterworth / Chebyshev design |

## Conventions used everywhere

- A sequence is a **(values, index) pair** `[x, n]` with contiguous, ascending indices.
- The math is **0-based**, MATLAB arrays are **1-based**, so every array access gets `+1`.
- Outputs carry their own index vector, so the time (or frequency) axis is always explicit.

## Files

**Module 1: sequences and convolution**

| File | Purpose |
|---|---|
| `signal_generator.m` | Interactive driver: asks for a signal type and parameters, plots a 3×3 `stem` grid |
| `gen_seq.m` | Sequence generator, returns values `x` and index vector `n` |
| `linear_conv.m` | Linear convolution from the definition, with index tracking |
| `circular_convolution.m` | N-point circular convolution, with index tracking |
| `test.m` | Compares both against `conv` / `cconv`, and shows linear convolution recovered through circular convolution |

**Module 2: DFT**

| File | Purpose |
|---|---|
| `dft.m` | M-point DFT from the definition (zero-pads if M > N, folds if M < N) |
| `idft.m` | M-point inverse DFT from the definition |
| `test_dft.m` | Experiments 1 to 3: DFT against `fft`, and IDFT recovery for M = N, M > N, M < N |
| `testdft2.m` | Experiment 5: folded DFT vs truncated DFT vs samples of the DTFT |

Result figures are in `results/`.

---

# Module 1: Sequences and convolution

## Signal generator

```matlab
[x, n] = gen_seq(type, n1, n2, varargin)     % n = n1:n2
```

| type | parameters (after `n1, n2`) | definition |
|---|---|---|
| `impulse`, `step` | `n0, A` | `A·δ(n−n0)`, `A·u(n−n0)` |
| `ramp` | `n0, A` | `A·(n−n0)·u(n−n0)` |
| `rect` | `n0, W, A` | `A` for `n0 ≤ n ≤ n0+W` (W+1 samples) |
| `sin` | `n0, A, N` | `A·sin(2π(n−n0)/N)` |
| `exp` | `n0, A, N` | `A·e^{j2π(n−n0)/N}` |
| `tri` | `n0, A, W` | `A·max(1 − \|n−n0\|/W, 0)` |
| `sinc` | `n0, A, N` | `A·sin(ωc m)/(πm)`, `m = n−n0`, `ωc = 2π/N`, value `A·ωc/π` at `m = 0` |
| `sign` | none | `+1` for `n ≥ 0`, `−1` for `n < 0` |

The `sinc` case avoids 0/0 at the centre by writing the limit value explicitly. Choosing a finite range `n1:n2` truncates the ideal low-pass response to an FIR filter.

## Linear convolution

`y(n) = Σₖ x(k) h(n−k)`. With 0-based positions, output position `s` only involves `k` from `max(0, s−Lh+1)` to `min(Lx−1, s)`. The real indices `nx`, `nh` are used only to label the output.

```matlab
Ly = Lx + Lh - 1;
ny = (min(nx)+min(nh)) : (max(nx)+max(nh));     % starts at nx(1)+nh(1)
y  = zeros(1, Ly);
for i = 1:Ly
    s = i - 1;
    for j = max(0, s-Lh+1) : min(Lx-1, s)
        y(i) = y(i) + x(j+1) * h(s-j+1);
    end
end
```

## Circular convolution

`y(n) = Σₖ x(k) h((n−k) mod N)`, with `N = max(Lx, Lh)`.

```matlab
N = max(Lx, Lh);
x = [x, zeros(1, N-Lx)];     % pad at the END only (front padding would shift the start index)
h = [h, zeros(1, N-Lh)];
ny = (min(nx)+min(nh)) : (min(nx)+min(nh)+N-1);
for i = 1:N
    s = i - 1;
    for j = 0:N-1
        y(i) = y(i) + x(j+1) * h(mod(s-j, N) + 1);
    end
end
```

- **Indexing rule.** If `x` starts at `a` and `h` at `b`, the 0-indexed result sits at indices `a+b … a+b+N−1`. The output is N-periodic, so this is one window of it.
- **Aliasing.** If `N < Lx+Lh−1`, the tail of the linear result wraps onto the front. If `N ≥ Lx+Lh−1` (pad both at the end), circular equals linear convolution with the same index vector.


| Result | Expected values | Index |
|---|---|---|
| Linear (8 points) | `2 6 5 5 6 3 2 1` | `-2 … 5` |
| Circular, N = 5 | `5 8 6 5 6` | `-2 … 2` |

Both outputs are compared against `conv` and `cconv`. As a check, each sums to 30 = sum(x)·sum(h), and the circular result is the linear result with its last three samples added onto the first three.

---

# Module 2: The DFT 

## Implementation

**DFT** of `M` points: `X(k) = Σₙ x(n)·W^(kn)`, with `W = e^(−j2π/M)`, `k = 0 … M−1`.

```matlab
% dft.m (core)
k = 0:M-1;
X = zeros(1, M);
W = exp(-1j*2*pi/M);
for i = 0:M-1
    for l = 0:M-1
        X(i+1) = X(i+1) + x(l+1) * W.^(i*l);
    end
end
```

**IDFT**: `x(n) = (1/M) Σₖ X(k)·W^(−kn)`, the same loop with the conjugate twiddle and a final `1/M`.

```matlab
% idft.m (core)
W = exp(1j*2*pi/M);
for i = 0:M-1
    for k = 0:M-1
        x(i+1) = x(i+1) + X(k+1) * W.^(k*i);
    end
    x(i+1) = x(i+1) / M;
end
```

**Choosing the DFT size `M` for an `N`-point input.** `dft(x, M)` handles all three cases:

| Case | What `dft.m` does |
|---|---|
| `M > N` | Zero-pads `x` at the end to length `M` |
| `M = N` | Uses `x` as it is |
| `M < N` | **Folds** the extra samples back onto the first ones (time-domain aliasing), so sample `n` and sample `n+M` are added |

```matlab
% folding for M < N
for i = 1:N-M
    x(i) = x(M+i) + x(i);
end
x = x(1:M);
```

Folding is what the DFT definition does automatically if the extra samples are left in. The twiddle `W^(kn)` repeats with period `M`, so samples `n` and `n+M` are indistinguishable. **Truncation** is a different choice: pass `x(1:M)` so the extra samples are discarded.

> MATLAB's `fft(x, M)` **truncates** when `M < N`. It will therefore match truncation, and will not match the folded `dft(x, M)`.

## Relationship to the DTFT

The DTFT of a finite sequence is `X(ω) = Σₙ x(n)·e^(−jωn)`. The M-point DFT is **M equally spaced samples of it**:

```
X[k] = X(ω) at ω = 2πk/M
```

Two consequences are tested in the experiments below: more zero-padding gives a denser sampling of the same curve, and folding (not truncation) keeps the samples on the DTFT of the *original* sequence.

## Test signal

```matlab
x = [1 2 0 1 5 6 1 8 5 9 12 10];     % N = 12, sum = 60
```

Every experiment uses this one signal, so the cases can be compared directly. In each, the DFT from `dft.m` is compared with `fft` (reference only), and the sequence recovered by `idft.m` is compared with the original.

---

## Experiments and results

### Case 1: M = N (12-point DFT)

```matlab
[X, k]   = dft(x, 12);
X2       = fft(x, 12);
[x_r, n] = idft(X);
```

![Case 1: M = N](results/Figure_2_DFT1.png)

- The real and imaginary parts of `dft` and `fft` coincide bin by bin (bin 0 equals the sum of the samples, 60).
- The IDFT returns all 12 samples and **recovers x(n) exactly**.

### Case 2: M > N (24-point DFT, zero-padded)

```matlab
[X, k]   = dft(x, 24);
X2       = fft(x, 24);
[x_r, n] = idft(X);
```

![Case 2: M > N](results/Figure_2_DFT2.png)

- `dft` and `fft` coincide on all 24 bins.
- The IDFT returns 24 samples: the original 12 followed by 12 zeros. It **recovers x(n)** plus the padding.
- Zero-padding adds no information. The extra bins are additional samples of the same DTFT curve.

### Case 3: M < N, truncated (8-point DFT of the first 8 samples)

```matlab
[X, k]   = dft(x(1:8), 8);      % extra samples discarded before the DFT
X2       = fft(x, 8);           % fft(x,8) truncates as well
[x_r, n] = idft(X);
```

![Case 3: M < N, truncation](results/Figure_2_DFT3.png)

- `dft` and `fft` coincide (bin 0 is 24, the sum of the first 8 samples, not 60).
- The IDFT recovers `[1 2 0 1 5 6 1 8]`, which is **x(n) for n = 0 to 7 only**. Samples 8 to 11 are lost.
- This spectrum belongs to a cropped signal, not to the original sequence.

### Case 4: M < N, folded (8-point DFT, extra samples added back)

```matlab
[X, k]   = dft(x, 8);           % folds samples 8..11 onto 0..3
[x_r, n] = idft(X);
```

Folding gives `x(0:3) + x(8:11) = [1+5, 2+9, 0+12, 1+10] = [6 11 12 11]`, so the sequence seen by the DFT is

```
[6 11 12 11 5 6 1 8]
```

- Bin 0 is 60, the sum of **all** 12 samples (truncation gave 24).
- The IDFT returns exactly the folded sequence, **not** the original x(n). Many different 12-sample sequences fold to the same 8 values, so the original cannot be recovered. The DFT is lossy whenever M < N.

### Case 5: Do the folded DFT samples lie on the original DTFT?

```matlab
[Xf, k] = dft(x, 8);                 % folded
[Xt, k] = dft(x(1:8), 8);            % truncated
wk = 2*pi*k/8;                       % DTFT of the ORIGINAL x at these frequencies
Xd = zeros(1, 8);
for n = 0:N-1
    Xd = Xd + x(n+1)*exp(-1j*wk*n);
end
```

The DTFT of the original sequence is also computed on a dense grid (1000 points over 0 to 2π) to draw the continuous curve.

![Case 5: folded vs truncated vs DTFT](results/Figure_3_truncated_vs_Folded_comparision.png)

- **Folded DFT (green) coincides with the DTFT samples (black)** in both the real and imaginary parts, and its magnitude lands on the DTFT curve at every bin.
- **Truncated DFT (red) does not.** It samples the DTFT of a different, cropped signal. At bin 0 the true value is 60, folding gives 60 and truncation gives 24. The missing 36 is exactly the sum of the discarded samples, 5+9+12+10.

To print the numeric difference, add these lines to `testdft2.m`:

```matlab
fprintf("max |folded    - DTFT| = %.2e\n", max(abs(Xf - Xd)));
fprintf("max |truncated - DTFT| = %.2e\n", max(abs(Xt - Xd)));
```

### Summary of the five cases

| Case | M | What the DFT uses | IDFT output | Recovers x(n)? | Spectrum samples |
|---|---|---|---|---|---|
| 1 | 12 (= N) | x | x | Yes | DTFT of x, 12 points |
| 2 | 24 (> N) | x + 12 zeros | x + 12 zeros | Yes (plus zeros) | DTFT of x, denser grid |
| 3 | 8 (< N) | x(0..7), tail discarded | x(0..7) | Only n = 0..7 | DTFT of the cropped signal |
| 4 | 8 (< N) | folded sequence | `[6 11 12 11 5 6 1 8]` | No | DTFT of x, 8 points |
| 5 | 8 (< N) | folded vs truncated | n/a | n/a | folded lies on the DTFT, truncated does not |

## Key findings

1. **The M-point DFT is M samples of the DTFT** at ω = 2πk/M. A different M samples the same curve at different frequencies, so plot against ω, not the bin index, when comparing sizes.
2. **M ≥ N loses nothing.** The IDFT returns x (plus zeros), and the extra bins from zero-padding only interpolate.
3. **M < N always loses information.** The two ways of reducing the length mean different things:
   - *Truncation* keeps exact samples of the original signal but analyses a cropped signal.
   - *Folding* keeps exact samples of the original signal's DTFT but mixes time regions hence ot invertible.
4. **`fft(x, M)` truncates for M < N**, so compare folded results against the DTFT, not against `fft`.

## limitations

- The folding loop in `dft.m` adds each extra sample once, which is correct when `N ≤ 2M`. For `N > 2M` a sample wraps more than once, and the general form is `x(mod(i-1, M) + 1) = x(mod(i-1, M) + 1) + x(i)` for `i = M+1:N`.
- Complexity is O(N²) complex multiplications, since the DFT is evaluated directly from the definition.
- The IDFT output is complex with negligible imaginary parts (rounding level) for a real signal. The plots show its real part.

