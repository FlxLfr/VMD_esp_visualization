# ==============================================================
# esp.tcl - ESP on the electron density isosurface
# Written by xyzToCubeToVMDVis.py (@@STAMP@@) from @@SOURCES@@
#
# Start:  vmd -e esp.tcl        or in the Tk Console:  source esp.tcl
# Commands: esp_view pi|edge|sigma, esp_iso, esp_range, esp_opacity, esp_snapshot
#
# This scene deliberately does NOT draw a colour bar. VMD has no legend object;
# one would have to put the bar into a molecule ID of its own as graphics and
# decouple it from the mouse - that is fiddly and never looks quite right.
# For the output images render_espVMD.py makes the scale with matplotlib, as a
# separate PNG in images/. The PyMOL pipeline works the same way.
# ==============================================================

# VMD versions know different display/material options. An unknown command
# raises an error, and an error aborts "source" then and there - everything
# after it is silently missing. Optional commands therefore go through _try.
proc _try {args} {
    if {[catch {uplevel 1 $args} err]} {
        puts "! skipped: $args   ($err)"
    }
}

# Colour scale with any number of anchor colours.
#
# "color scale method" knows only three: lowest value, middle, highest. The
# rainbow ramp of the PyMOL pipeline has five (red, yellow, green, cyan,
# blue), and with three of them it turns into red-green-blue with olive and
# teal in between - a different ramp, not the same one.
#
# VMD's colour scale, however, is a table of 1024 colours whose IDs sit behind
# the named colours ([colorinfo num], normally 33). Each one can be set
# individually with "color change rgb". That makes the ramp freely choosable
# and exactly the same as over there.
#
# Important: call this AFTER "color scale method". A later method command
# rebuilds the table and overwrites everything set here.
proc esp_ramp {stops} {
    if {[llength $stops] < 2} { return }        ;# empty = VMD's own ramp
    set base [colorinfo num]
    set n    1024
    set segs [expr {[llength $stops] - 1}]
    _try display update off
    for {set i 0} {$i < $n} {incr i} {
        set t [expr {double($i) / ($n - 1) * $segs}]
        set k [expr {int($t)}]
        if {$k >= $segs} { set k [expr {$segs - 1}] }
        set f  [expr {$t - $k}]
        set c0 [lindex $stops $k]
        set c1 [lindex $stops [expr {$k + 1}]]
        color change rgb [expr {$base + $i}] \
            [expr {[lindex $c0 0] + $f * ([lindex $c1 0] - [lindex $c0 0])}] \
            [expr {[lindex $c0 1] + $f * ([lindex $c1 1] - [lindex $c0 1])}] \
            [expr {[lindex $c0 2] + $f * ([lindex $c1 2] - [lindex $c0 2])}]
    }
    _try display update on
}

# --- 1) Data --------------------------------------------------
# The cube file brings the atoms with it, in Bohr and from the same source as
# the grid - so a separate structure file cannot slip out of register.
mol delete all
set espmol [mol new @@RHO_CUBE@@ type cube waitfor all]
mol addfile @@ESP_CUBE@@ type cube waitfor all

# The second cube file appends an identical set of coordinates as a frame.
if {[molinfo $espmol get numframes] > 1} {
    animate delete beg 1 end [expr {[molinfo $espmol get numframes] - 1}] $espmol
}

set VOL_RHO @@VOL_RHO@@ ;# electron density
set VOL_ESP @@VOL_ESP@@ ;# electrostatic potential
set ISO     @@ISO@@
set RANGE   @@RANGE@@
set SCALE   @@SCALE@@   ;# a number or "auto"
set FILL    @@FILL@@    ;# fraction of the window height the molecule should fill

# --- 2) Representations ---------------------------------------
# THE ORDER IS NOT COSMETIC: VMD draws reps by number and writes to the depth
# buffer for transparent surfaces too. With the isosurface first, the skeleton
# drawn afterwards fails the depth test everywhere - a visibly translucent
# surface with nothing to be seen behind it.
# So: the opaque skeleton first, then the transparent surface.
mol delrep 0 $espmol
set REP_MOL  0
set REP_SURF 1

_try color Name C gray
# Licorice <bond radius> <sphere resolution> <cylinder resolution>. The two
# resolutions are facet counts, not sizes: 24 keeps the cylinders round at
# full zoom without costing noticeable render time at this image size.
mol representation Licorice @@STICK@@ 24.000000 24.000000
mol color Name
mol selection {all}
mol material Opaque
mol addrep $espmol

# Isosurface rho = ISO, coloured by the SECOND volumetric dataset. This is the
# VMD equivalent of ramp_new + set surface_color in PyMOL: VMD has no ramp
# object, it colours directly by the values of another grid in the same
# molecule ID.
#   Isosurface <iso> <volID> <show 0=surface> <draw 0=solid> <step> <size>
mol representation Isosurface $ISO $VOL_RHO 0 0 1 1
mol color Volume $VOL_ESP
mol selection {all}
mol material Transparent
mol addrep $espmol

# @@RANGE_NEG@@ .. +@@RANGE@@ a.u. = @@KJ@@ kJ/(mol*e). RWB: red negative.
color scale method @@COLORSCALE@@
color scale midpoint 0.5
mol scaleminmax $espmol $REP_SURF @@RANGE_NEG@@ @@RANGE@@

# Anchor colours of the ramp, taken over from the PyMOL pipeline. Empty means:
# VMD's built-in scale @@COLORSCALE@@ stays as it is.
set RAMP_STOPS {@@RAMP_STOPS@@}
esp_ramp $RAMP_STOPS

# --- 3) Commands ----------------------------------------------

# Zoom from the size of the molecule instead of from a guessed number.
#
# The visible part of the world is about "display height"/2 tall - measured
# empirically against a real Tachyon render, VMD does not report the size
# anywhere directly. For the molecule to fill FILL of the image height it must
# therefore hold that:  scale = FILL * (height/2) / (2*r).
#
# r is the radius of the bounding sphere plus about 1.8 A for the isosurface -
# via the sphere, so that none of the three views is cut off.
proc esp_fit {} {
    global espmol SCALE FILL
    if {$SCALE ne "auto"} { scale to $SCALE ; return }
    set sel [atomselect $espmol all]
    set mm [measure minmax $sel]
    $sel delete
    set ex [expr {[lindex [lindex $mm 1] 0] - [lindex [lindex $mm 0] 0]}]
    set ey [expr {[lindex [lindex $mm 1] 1] - [lindex [lindex $mm 0] 1]}]
    set ez [expr {[lindex [lindex $mm 1] 2] - [lindex [lindex $mm 0] 2]}]
    set r  [expr {0.5 * sqrt($ex*$ex + $ey*$ey + $ez*$ez) + 1.8}]
    # Fallback in case this VMD version does not hand out the height: 6.0 is
    # the default. Without that guard an empty return value takes the whole
    # scene down with it, because esp_fit is already called at start-up.
    if {[catch {set h [display get height]}] || ![string is double -strict $h]
        || $h <= 0} {
        set h 6.0
    }
    scale to [expr {$FILL * 0.25 * $h / $r}]
}

# Centre on the ATOMS, not on the grid box: resetview fits the view to all
# reps, and the isosurface carries the extent of the whole grid.
proc esp_center {} {
    global espmol
    set sel [atomselect $espmol all]
    set c [measure center $sel]
    molinfo $espmol set center_matrix [list [transoffset [vecscale -1.0 $c]]]
    molinfo $espmol set rotate_matrix [list [transidentity]]
    molinfo $espmol set global_matrix [list [transidentity]]
    molinfo $espmol set scale_matrix  [list [transidentity]]
    $sel delete
    esp_fit
}

# The three views sit here as finished rotation matrices.
# xyzToCubeToVMDVis.py computed them from the geometry: the ring normal for
# pi, the true carbon-halogen axis for sigma - the same calculation as in the
# PyMOL pipeline, so that both image sets really show the same thing.
#
# This used to read "pi, then rotate x by -90". That assumed the molecule lies
# planar in the xy plane with the C-X axis pointing along -y. For the
# halobenzenes from Turbomole that held; for the substituted pyridines it did
# not, and there sigma looked past the hole. It would have been almost
# impossible to notice: in every direction there is a round coloured surface,
# and the far side looks like a sigma hole.
#
# sigma axis of this molecule: @@AXIS_LABEL@@
set ROT_PI    @@ROT_PI@@
set ROT_EDGE  @@ROT_EDGE@@
set ROT_SIGMA @@ROT_SIGMA@@

proc esp_view {which} {
    global espmol ROT_PI ROT_EDGE ROT_SIGMA
    esp_center
    switch -- $which {
        pi      { set m $ROT_PI }
        edge    { set m $ROT_EDGE }
        sigma   { set m $ROT_SIGMA }
        default { puts "esp_view: pi | edge | sigma" ; return }
    }
    molinfo $espmol set rotate_matrix [list $m]
    if {$which eq "sigma"} {
        puts "Note: in the axial view the transparent layers"
        puts "  stack up - esp_opacity 1.0 or esp_snapshot"
        puts "  (Tachyon) gives a clean image."
    }
    puts "View: $which"
}

proc esp_iso {value} {
    global espmol VOL_RHO REP_SURF ISO
    set ISO $value
    mol modstyle $REP_SURF $espmol Isosurface $value $VOL_RHO 0 0 1 1
    puts "Isovalue: rho = $value a.u."
}

proc esp_range {half} {
    global espmol REP_SURF RANGE
    set RANGE $half
    mol scaleminmax $espmol $REP_SURF [expr {-1.0 * $half}] $half
    puts "Colour scale: +/- $half a.u."
}

# The line of sight passes through TWO layers of the closed surface; at opacity
# a, (1-a)^2 survives. 0.5 is the workable compromise, 1.0 shows only the
# surface, below 0.3 the colour becomes too pale.
proc esp_opacity {value} {
    _try material change opacity Transparent $value
    puts "Opacity: $value"
}

# Ray traced. Fast alternative: render snapshot $name.png
proc esp_snapshot {name} {
    render TachyonInternal $name.tga
    puts "written: $name.tga"
}

# --- 4) Display -----------------------------------------------
_try display projection Orthographic
_try display depthcue off
_try display shadows off
_try display culling off
_try axes location off
_try color Display Background white
# Without this the near clipping plane cuts a slice out when zooming.
_try display nearclip set 0.010000
# True transparency instead of a dither pattern ("screen door"); old drivers
# cannot do GLSL - then the pattern stays, the Tachyon images are clean anyway.
_try display rendermode GLSL
_try material change opacity  Transparent @@OPACITY@@
_try material change diffuse  Transparent 0.750000
# Keep the specular highlight small - a milky veil eats up the view through.
_try material change specular Transparent 0.100000
_try material change shininess Transparent 0.300000

# --- 5) Start -------------------------------------------------
esp_view pi

puts ""
puts "esp.tcl loaded. rho = $ISO a.u., colour scale +/- $RANGE a.u.@@STATS@@"
puts ""
