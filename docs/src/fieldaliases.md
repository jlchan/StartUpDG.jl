# Fieldname aliases

`Base.getproperty` is overloaded for both `RefElemData` and `MeshData`, so you can reference `rd.Dr` instead of `md.Drst[1]` or `md.rxJ` instead of `md.rstxyzJ[1,1]`. `@unpack` also works with aliased fieldnames.

The aliases for `RefElemData` are below:
1. `r` for `rst[1]` (similarly for `s` and `t`)
2. `rq` for `rstq[1]` (similarly for `sq` and `tq`)
3. `rf` for `rstf[1]` (similarly for `sf` and `tf`)
4. `rp` for `rstp[1]` (similarly for `sp` and `tp`)
5. `Dr` for `Drst[1]` (similarly for `Ds` and `Dt`)
6. `Nfaces, Np, Nq, Nfq` for number of faces, volume nodes, quadrature points, and face quadrature points. 

The aliases for `MeshData` are below:
1. `VX` for `VXYZ[1]` (similarly for `VY`, `VZ`)
2. `x` for `xyz[1]` (similarly for `y`, `z`)
3. `xf` for `xyzf[1]` (similarly for `yf`, `zf`)
4. `xq` for `xyzq[1]` (similarly for `yq`, `zq`)
5. `nxJ` for `nxyzJ[1]` (similarly for `nyJ`, `nzJ`)
6. `rxJ, sxJ` for `rstxyzJ[1,1],rstxyzJ[1,2]` (similarly for `txJ,ryJ,syJ,tyJ,rzJ,szJ,tzJ`)
