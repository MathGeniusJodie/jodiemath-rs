/// Every public `fn(f32) -> f32` in the crate, for the whole-API sweeps
/// (special_matrix, nan_payload, worst_corpus). `assert_unary_fns_complete`
/// fails loudly when src/lib.rs gains or loses one and this list has not
/// followed, which is how the three copies this replaced drifted: 8 public
/// functions were missing and 8 deleted ones were still listed.
const UNARY_FNS: &[(&str, fn(f32) -> f32)] = &[
    ("acos", acos),
    ("acosd", acosd),
    ("acosh", acosh),
    ("acospi", acospi),
    ("asin", asin),
    ("asind", asind),
    ("asinh", asinh),
    ("asinpi", asinpi),
    ("atan", atan),
    ("atan_bounded", atan_bounded),
    ("atan_latency", atan_latency),
    ("atand", atand),
    ("atanh", atanh),
    ("cbrt", cbrt),
    ("cbrt_accurate", cbrt_accurate),
    ("cbrt_accurate_unchecked", cbrt_accurate_unchecked),
    ("cbrt_approx", cbrt_approx),
    ("cbrt_fast", cbrt_fast),
    ("cbrt_unchecked", cbrt_unchecked),
    ("cos", cos),
    ("cos_wide", cos_wide),
    ("cosh", cosh),
    ("cosh_checked", cosh_checked),
    ("cosh_narrow", cosh_narrow),
    ("cosh_throughput", cosh_throughput),
    ("coshm1", coshm1),
    ("cospi", cospi),
    ("erf", erf),
    ("erfc", erfc),
    ("erfc_inv", erfc_inv),
    ("erfinv", erfinv),
    ("exp", exp),
    ("exp10", exp10),
    ("exp10_checked", exp10_checked),
    ("exp10m1", exp10m1),
    ("exp2", exp2),
    ("exp2_approx", exp2_approx),
    ("exp2_checked", exp2_checked),
    ("exp2m1", exp2m1),
    ("exp_checked", exp_checked),
    ("exp_m1_over_x_narrow", exp_m1_over_x_narrow),
    ("exp_narrow", exp_narrow),
    ("expm1", expm1),
    ("expm1_checked", expm1_checked),
    ("expm1_narrow", expm1_narrow),
    ("fast_round_int", fast_round_int),
    ("gelu", gelu),
    ("linear_to_srgb", linear_to_srgb),
    ("ln", ln),
    ("ln_unchecked", ln_unchecked),
    ("log10", log10),
    ("log10_unchecked", log10_unchecked),
    ("log10p1", log10p1),
    ("log1p", log1p),
    ("log1pmx", log1pmx),
    ("log2_approx", log2_approx),
    ("log2p1", log2p1),
    ("log_2", log_2),
    ("log_2_unchecked", log_2_unchecked),
    ("logit", logit),
    ("logsigmoid", logsigmoid),
    ("logsigmoid_checked", logsigmoid_checked),
    ("pow_2_3", pow_2_3),
    ("pow_3_2", pow_3_2),
    ("rcbrt", rcbrt),
    ("rcp_approx", rcp_approx),
    ("rsqrt", rsqrt),
    ("rsqrt_approx", rsqrt_approx),
    ("sigmoid", sigmoid),
    ("sigmoid_fast", sigmoid_fast),
    ("sigmoid_grad", sigmoid_grad),
    ("silu", silu),
    ("silu_checked", silu_checked),
    ("sin", sin),
    ("sin_wide", sin_wide),
    ("sinc_unnormalized", sinc_unnormalized),
    ("sinh", sinh),
    ("sinh_checked", sinh_checked),
    ("sinh_narrow", sinh_narrow),
    ("sinh_throughput", sinh_throughput),
    ("sinpi", sinpi),
    ("softsign", softsign),
    ("sqrt1pm1", sqrt1pm1),
    ("sqrt_approx", sqrt_approx),
    ("srgb_to_linear", srgb_to_linear),
    ("tan", tan),
    ("tan_wide", tan_wide),
    ("tanh", tanh),
    ("tanh_grad", tanh_grad),
    ("tanpi", tanpi),
];

fn assert_unary_fns_complete() {
    let src = include_str!("../../src/lib.rs");
    let public: Vec<&str> = src
        .lines()
        .filter_map(|l| l.strip_prefix("pub fn "))
        .filter_map(|l| {
            let (name, rest) = l.split_once('(')?;
            let (arg, _) = rest.split_once(')')?;
            (arg.ends_with(": f32") && !arg.contains(',') && rest.contains(") -> f32 {"))
                .then_some(name)
        })
        .collect();
    let listed: Vec<&str> = UNARY_FNS.iter().map(|(n, _)| *n).collect();
    let missing: Vec<&&str> = public.iter().filter(|n| !listed.contains(n)).collect();
    let stale: Vec<&&str> = listed.iter().filter(|n| !public.contains(n)).collect();
    assert!(
        missing.is_empty() && stale.is_empty(),
        "examples/support/unary_fns.rs is out of date with src/lib.rs: missing {missing:?}, no longer public {stale:?}"
    );
}
