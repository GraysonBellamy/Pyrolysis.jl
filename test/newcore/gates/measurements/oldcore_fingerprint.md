# The old core's bit-level fingerprint

Produced by `test/newcore/gates/oldcore_fingerprint.jl` at commit 603e863. Julia 1.13.0. The hash of 3072 reaction rates and heats through `Physics.evaluate_reactions!` (orders 0, 0.5, 1, 2; depletion thresholds 1e-3, 0.1 and 10 kg/m³, the limiter on and off; 16 solid concentrations from −1 to 1e3 kg/m³ at 300, 450, 700 and 1200 K) and `hash((T, ξ, t))` of eight old-core solves at their CI resolutions, against the same digests on 17bc477, the tree before the finite-volume core's old-core edits (Julia 1.13.0).

| item | 17bc477 | this tree | equal | solve |
|---|---|---|---|---|
| reaction rates and heats (3072 values) | `0xc8fea394cd668f5a` | `0xc8fea394cd668f5a` | true |  |
| V-RX-01 | `0x029fd72814cbb78d` | `0x029fd72814cbb78d` | true | Success, 7 saved times |
| V-RX-02a | `0x1c0b673ed2de886a` | `0x1c0b673ed2de886a` | true | Success, 6 saved times |
| V-RX-02c | `0xd3d03eade4a8bd42` | `0xd3d03eade4a8bd42` | true | Success, 7 saved times |
| V-RX-02d | `0x552972dd23e92d46` | `0x552972dd23e92d46` | true | Success, 21 saved times |
| V-RX-03 | `0xdcb652b0a1de03d5` | `0xdcb652b0a1de03d5` | true | Success, 7 saved times |
| V-RX-04a | `0xba656d55a93f8d82` | `0xba656d55a93f8d82` | true | Success, 6 saved times |
| V-CPL-01 | `0x98d675748199bed5` | `0x98d675748199bed5` | true | Success, 31 saved times |
| V-CPL-04 | `0x2765662b2c3e4833` | `0x2765662b2c3e4833` | true | Success, 11 saved times |

Every digest equals its reference: the old core is bit-identical.
