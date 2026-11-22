/* BCHVision engine: GF(2^m), binary narrow-sense primitive BCH, BM + Chien decoder. Polys are arrays, index = power of x. */
(function (root) {
  const PP = { 3: 11, 4: 19, 5: 37, 6: 67, 7: 137, 8: 285 };
  function field(m) {
    const n = (1 << m) - 1, exp = new Array(2 * n), log = new Array(n + 1);
    let x = 1;
    for (let i = 0; i < n; i++) { exp[i] = exp[i + n] = x; log[x] = i; x <<= 1; if (x >> m) x ^= PP[m]; }
    return { m, n, exp, log, pp: PP[m],
      mul: (a, b) => (a && b ? exp[log[a] + log[b]] : 0),
      inv: a => exp[(n - log[a]) % n] };
  }
  function cosets(n, t) {
    const seen = new Set(), out = [];
    for (let i = 1; i <= 2 * t; i++) {
      if (i >= n || seen.has(i)) continue;
      const c = []; let j = i;
      do { c.push(j); seen.add(j); j = (j * 2) % n; } while (j !== i);
      out.push(c);
    }
    return out;
  }
  function capability(m) {
    const n = (1 << m) - 1, rows = [];
    for (let t = 1; 2 * t < n; t++) {
      const deg = cosets(n, t).reduce((s, c) => s + c.length, 0);
      if (n - deg <= 0) break;
      rows.push({ t, n, k: n - deg, deg });
    }
    return rows;
  }
  function build(m, t) {
    const F = field(m), n = F.n, cs = cosets(n, t);
    const pmul = (a, b) => { const r = new Array(a.length + b.length - 1).fill(0);
      a.forEach((x, i) => b.forEach((y, j) => { r[i + j] ^= F.mul(x, y); })); return r; };
    const minpolys = cs.map(c => c.reduce((p, j) => pmul(p, [F.exp[j], 1]), [1]));
    const g = minpolys.reduce((p, q) => pmul(p, q), [1]);
    return { F, n, t, m, cosets: cs, minpolys, g, deg: g.length - 1, k: n - (g.length - 1), delta: 2 * t + 1 };
  }
  function encode(C, msg) { // msg: k bits, index = power of x
    const { n, k, deg, g } = C, c = new Array(n).fill(0);
    for (let i = 0; i < k; i++) c[deg + i] = msg[i];
    const r = c.slice();
    for (let j = n - 1; j >= deg; j--) if (r[j]) for (let i = 0; i <= deg; i++) r[j - deg + i] ^= g[i];
    const rem = r.slice(0, deg);
    for (let i = 0; i < deg; i++) c[i] = rem[i];
    return { cw: c, rem };
  }
  function syndromes(C, r) {
    const S = [];
    for (let i = 1; i <= 2 * C.t; i++) { let s = 0; for (let j = 0; j < C.n; j++) if (r[j]) s ^= C.F.exp[(i * j) % C.n]; S.push(s); }
    return S;
  }
  function decode(C, r) {
    const { F, n, t } = C, { exp, mul, inv } = F, S = syndromes(C, r);
    if (S.every(s => !s)) return { S, lam: [1], pos: [], ok: true, cw: r.slice() };
    let Cc = [1], B = [1], L = 0, mm = 1, b = 1;
    for (let k = 0; k < 2 * t; k++) {
      let d = S[k]; for (let i = 1; i <= L; i++) d ^= mul(Cc[i] || 0, S[k - i]);
      if (!d) { mm++; continue; }
      const T = Cc.slice(), co = mul(d, inv(b));
      for (let i = 0; i < B.length; i++) Cc[i + mm] = (Cc[i + mm] || 0) ^ mul(co, B[i]);
      if (2 * L <= k) { L = k + 1 - L; B = T; b = d; mm = 1; } else mm++;
    }
    const lam = Array.from({ length: L + 1 }, (_, i) => Cc[i] || 0), pos = [];
    for (let j = 0; j < n; j++) { let v = 0; lam.forEach((c, i) => { if (c) v ^= mul(c, exp[((n - j) * i) % n]); }); if (!v) pos.push(j); }
    const cw = r.slice(); pos.forEach(j => { cw[j] ^= 1; });
    const ok = L <= t && pos.length === L && syndromes(C, cw).every(s => !s);
    return { S, lam, pos, ok, cw: ok ? cw : r.slice() };
  }
  function dmin(C) { // exhaustive, only when small
    const { n, k, deg, g } = C;
    if (n > 31 || k > 20) return null;
    const rows = []; for (let i = 0; i < k; i++) { let v = 0; g.forEach((b, p) => { if (b) v |= 1 << (p + i); }); rows.push(v >>> 0); }
    const pc = v => { v -= (v >>> 1) & 0x55555555; v = (v & 0x33333333) + ((v >>> 2) & 0x33333333); return (((v + (v >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24; };
    let cur = 0, best = n;
    for (let i = 1; i < 2 ** k; i++) { let z = 0, x = i; while (!(x & 1)) { x >>= 1; z++; } cur = (cur ^ rows[z]) >>> 0; const w = pc(cur); if (w < best) best = w; }
    return best;
  }
  const api = { field, cosets, capability, build, encode, decode, syndromes, dmin, PP };
  if (typeof module !== 'undefined') module.exports = api; else root.BCH = api;
})(this);