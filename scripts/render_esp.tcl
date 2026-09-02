# ==============================================================
# render_esp.tcl - the three standard views, ray traced, into images/
#
# Call it from the molecule folder (where esp.tcl and the cubes live):
#
#     vmd -e ../../scripts/render_esp.tcl
#
# Normally render_espVMD.py calls this script and takes care of the colour
# bar and settings.txt afterwards. Started directly it only does the VMD part.
#
# Can be set beforehand (in the Tk Console, then "source"):
#   set ESP_RES    {1600 1280}   window size = image size
#   set ESP_QUIT   0             leave VMD open after rendering
#   set ESP_OPAQUE 1             render the isosurface opaque
#   set ESP_BG     black          background colour
#   set ESP_VIEWS  {sigma}       single views only
#   set ESP_SNAPSHOT 1           window capture instead of ray tracing
#   set ESP_SCENE  esp_check.tcl a different scene file (self test)
#   set ESP_OUTDIR images_check   a different target folder (self test)
#   set ESP_PREFIX brombenzol     file prefix (default: folder name)
# ==============================================================

if {![info exists ESP_RES]}    { set ESP_RES {1600 1280} }
if {![info exists ESP_QUIT]}   { set ESP_QUIT 1 }
if {![info exists ESP_OPAQUE]}  { set ESP_OPAQUE 0 }
if {![info exists ESP_BG]}     { set ESP_BG white }
# Appended to the file name - with several background colours in one run it
# tells the image sets apart (<molecule>_pi_black.tga).
if {![info exists ESP_SUFFIX]} { set ESP_SUFFIX "" }
if {![info exists ESP_VIEWS]}  { set ESP_VIEWS {pi edge sigma} }
if {![info exists ESP_SNAPSHOT]} { set ESP_SNAPSHOT 0 }
# The self test writes its scene as esp_check.tcl, so that it does not
# overwrite the committed esp.tcl of a real run.
if {![info exists ESP_SCENE]}  { set ESP_SCENE esp.tcl }
# Target folder and prefix come from render_espVMD.py. If they are not set,
# the same defaults apply as for a call by hand: images/ and the name of the
# molecule folder.
if {![info exists ESP_OUTDIR]} { set ESP_OUTDIR images }
if {![info exists ESP_PREFIX]} { set ESP_PREFIX [file tail [pwd]] }

if {![file exists $ESP_SCENE]} {
    puts "render_esp.tcl: $ESP_SCENE not found."
    puts "  Start it from the molecule folder, or run xyzToCubeToVMDVis.py first."
    if {$ESP_QUIT} { quit }
    return
}

source $ESP_SCENE

# Emergency exit for the axial view: there every line of sight crosses many
# transparent layers of the isosurface, and together with ambient occlusion
# that crashes Tachyon. Rendered opaque, the transparency recursion goes away.
# render_espVMD.py only switches this on when it was needed.
if {$ESP_OPAQUE} { esp_opacity 1.0 }

# --- Render quality -------------------------------------------
# Drop shadows stay off, without an option. They cast the licorice sticks onto
# the isosurface as grey capsules: in the pi image an offset double sits behind
# every stick, in the axial view a stick shadow lies in the middle of the
# sphere. Both look like an artefact of the data and are none, and the PyMOL
# pipeline sets ray_shadows 0 for the same reason - which is what keeps the two
# image sets comparable.
#
# Ambient occlusion is off as well, and there is no switch for it either. It
# darkens the hollows between the CH bulges, which does look richer, but it
# also lines the sticks with grey doubles on the isosurface, it has no
# counterpart in the PyMOL pipeline, and combined with many transparent layers
# it was the most reliable way to make Tachyon abort in the axial view.
# Anyone who wants to see either effect switches it on in the Tk Console:
# 'display shadows on' or 'display ambientocclusion on', then
# 'esp_snapshot <name>'.
_try display shadows off
_try display ambientocclusion off
_try color Display Background $ESP_BG
_try display antialias on
_try display depthcue off
_try display resize [lindex $ESP_RES 0] [lindex $ESP_RES 1]
display update

# TachyonInternal renders at window size. If the screen is smaller than
# ESP_RES, Windows clamps the window and the image comes out correspondingly
# smaller - which is why the actual size is reported.
puts "Window: [display get size]" ; flush stdout

# The window capture is the emergency exit for the axial view: Tachyon bails
# out there on the number of transparent layers crossed, while OpenGL draws
# them without recursion. The image is a little less fine but keeps the
# transparency - better than an opaque substitute.
# Careful: the capture copies the screen contents, so the VMD window must not
# be covered.
#
# The renderer is hard-wired, and the fastest one available is deliberately
# not the one used. There is no option for it - see below.
#
# VMD does not colour a "mol color Volume" isosurface vertex by vertex, it
# lays the 1024 entries of the colour scale onto the surface as a 1D TEXTURE.
# A renderer that does not export textures therefore draws the surface in the
# plain material colour - white, in this scene. TachyonLOptiXInternal is such
# a renderer. VMD offers it on every machine with an NVIDIA card, it is a few
# tenths of a second faster over three images, and it says so itself while it
# happens:
#
#     Warning) Texture mapping not exported for this renderer
#     Warning) Unimplemented features may negatively affect the appearance
#
# The images are then geometrically correct and completely colourless - the
# one thing the whole pipeline is about is missing. This used to be selected
# automatically from what the machine offered, which made the picture depend
# on the graphics card of whoever ran the pipeline; on an NVIDIA box the whole
# image set came out white. TachyonInternal is in every VMD build, needs no
# GPU and takes about a second per image here, so there is nothing to gain
# from making this selectable - a second renderer could only reintroduce the
# same bug.
set RENDERER TachyonInternal
if {$ESP_SNAPSHOT} { set RENDERER snapshot }
puts "Renderer: $RENDERER" ; flush stdout

# --- Render ---------------------------------------------------
file mkdir $ESP_OUTDIR
set prefix $ESP_PREFIX
puts "Target: $ESP_OUTDIR" ; flush stdout

# Guard every view separately: if one aborts - Tachyon can bail out on large
# isosurfaces with ambient occlusion - the others should still be produced.
# And flush after every line, so that the log shows how far VMD got if the
# process dies halfway through.
foreach view $ESP_VIEWS {
    set out [file join $ESP_OUTDIR "${prefix}_${view}${ESP_SUFFIX}.tga"]
    puts "== $view: start ==" ; flush stdout
    if {[catch {esp_view $view ; display update} err]} {
        puts "! $view: view failed: $err" ; flush stdout
        continue
    }
    if {[catch {render $RENDERER $out} err]} {
        puts "! $view: render failed: $err" ; flush stdout
    } else {
        puts "-> $out" ; flush stdout
    }
}

puts "render_esp.tcl done. render_espVMD.py makes the PNGs and the colour bar."
if {$ESP_QUIT} { quit }
