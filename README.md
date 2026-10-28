# Discrete-Time Signal Processing from Scratch (MATLAB)

Design-level implementation of core DSP building blocks **without** calling `conv`, `cconv`, `fft` or `filter`. Built-in functions are used only as reference oracles in the tests.

**Status:** Phase 1 (sequences, linear and circular convolution) is done. Phase 2 (DFT/FFT, DTFT, windowed sinc and filter design) is next.

## Files

| File | Purpose |
|---|---|
| `signal_generator.m` | Interactive driver: asks for signal type and parameters, plots in a 3×3 `stem` grid |
| `gen_seq.m` | Sequence generator: returns values `x` and index vector `n` |
| `linear_conv.m` | Linear convolution from the definition, with index tracking |
| `circular_convolution.m` | N-point circular convolution, with index tracking |
| `test.m` | Compares both against `conv` / `cconv` and shows linear conv recovered from circular conv |

## Core convention

Every sequence is a **(values, index) pair**, `[x, n]`, with contiguous ascending indices. Math is 0-based and MATLAB arrays are 1-based, so every array access gets `+1`.

## 1. Signal generator (`gen_seq`)

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

The `sinc` case avoids the 0/0 at the centre by writing the limit value explicitly. Choosing a finite range `n1:n2` is what truncates the ideal low-pass response to an FIR filter.

```matlab
m  = n - n0;
x  = A * (wc/pi) * ones(size(n));          % limit value at m = 0
nz = (m ~= 0);
x(nz) = A .* sin(wc * m(nz)) ./ (pi * m(nz));
```

## 2. Linear convolution

y(n) = Σₖ x(k) h(n−k). With 0-based array positions, output position `s` is valid only when both `k` and `s−k` fall inside their arrays, so `k` runs from `max(0, s−Lh+1)` to `min(Lx−1, s)`. The real indices `nx`, `nh` are used only for the output label.

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

Output length is `Lx+Lh−1`, starting at `nx(1)+nh(1)`.

## 3. Circular convolution

y(n) = Σₖ x(k) h((n−k) mod N), for `N = max(Lx, Lh)`.

```matlab
N = max(Lx, Lh);
x = [x, zeros(1, N-Lx)];     % pad at the END only (front padding would change the start index)
h = [h, zeros(1, N-Lh)];
ny = (min(nx)+min(nh)) : (min(nx)+min(nh)+N-1);
for i = 1:N
    s = i - 1;
    for j = 0:N-1
        y(i) = y(i) + x(j+1) * h(mod(s-j, N) + 1);
    end
end
```

**Indexing rule.** If `x` starts at `a` and `h` at `b`, the 0-indexed array result `y0` sits at indices `a+b … a+b+N−1`. Because the output is N-periodic, this is one window of the periodic result. The same sequence on `0 … N−1` is `circshift(y0, a+b)`.

**Aliasing.** If `N < Lx+Lh−1`, the tail of the linear result wraps onto the front:

`y_circ(a+b+j) = Σₖ y_lin(a+b+j+kN)`

If `N ≥ Lx+Lh−1` (zero-pad both sequences at the end), circular convolution equals linear convolution with the same index vector.

## 4. Verification example

```matlab
x = [1 2 0 1 1];  nx = -2:2;
h = [2 2 1 1];    nh =  0:3;
```

| Result | Expected values | Index |
|---|---|---|
| Linear (8 points) | `2 6 5 5 6 3 2 1` | `-2 … 5` |
| Circular, N = 5 | `5 8 6 5 6` | `-2 … 2` |

Checks: both sums equal 30 = sum(x)·sum(h); the circular result is the linear result with its last three samples added onto the first three; `circular_convolution` on inputs zero-padded to length 8 equals `linear_conv`. Compare against `conv` and `cconv(x, h, 5)`.

## Next

- DFT from the definition, then the DFT matrix (its link to the circulant matrix explains why the DFT turns circular convolution into multiplication)
- Radix-2 FFT, Bluestein for arbitrary N, overlap-add convolution
- Windowed-sinc FIR design and DFT analysis of the ideal and truncated low-pass filters
- Butterworth and Chebyshev design with a quantization study
