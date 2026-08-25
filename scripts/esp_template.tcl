# ==============================================================
# esp.tcl - ESP auf der Elektronendichte-Isoflaeche
# Erzeugt von xyzToCubeToVMDVis.py (@@STAMP@@) aus @@SOURCES@@
#
# Start:  vmd -e esp.tcl        oder in der Tk Console:  source esp.tcl
# Befehle: esp_view pi|edge|sigma, esp_iso, esp_range, esp_opacity,
#          esp_colorbar, esp_colorbar_off, esp_snapshot
# ==============================================================

# VMD-Versionen kennen unterschiedliche display-/material-Optionen. Ein
# unbekannter Befehl wirft einen Fehler, und ein Fehler bricht "source" an Ort
# und Stelle ab - alles danach fehlt dann kommentarlos. Optionales laeuft
# deshalb ueber _try.
proc _try {args} {
    if {[catch {uplevel 1 $args} err]} {
        puts "! uebersprungen: $args   ($err)"
    }
}

# --- 1) Daten -------------------------------------------------
# Die Cube-Datei bringt die Atome selbst mit, in Bohr und aus derselben Quelle
# wie das Gitter - eine separate Strukturdatei kann also nicht verrutschen.
mol delete all
set espmol [mol new @@RHO_CUBE@@ type cube waitfor all]
mol addfile @@ESP_CUBE@@ type cube waitfor all

# Die zweite Cube-Datei haengt einen identischen Koordinatensatz als Frame an.
if {[molinfo $espmol get numframes] > 1} {
    animate delete beg 1 end [expr {[molinfo $espmol get numframes] - 1}] $espmol
}

set VOL_RHO @@VOL_RHO@@ ;# Elektronendichte
set VOL_ESP @@VOL_ESP@@ ;# elektrostatisches Potential
set ISO     @@ISO@@
set RANGE   @@RANGE@@
set SCALE   @@SCALE@@   ;# Zahl oder "auto"
set FILL    @@FILL@@    ;# Anteil der Fensterhoehe, den das Molekuel fuellen soll

# --- 2) Darstellungen -----------------------------------------
# REIHENFOLGE IST NICHT KOSMETIK: VMD zeichnet Reps nach Nummer und schreibt
# auch fuer transparente Flaechen in den Tiefenpuffer. Steht die Isoflaeche
# zuerst, faellt das danach gezeichnete Geruest ueberall aus dem Tiefentest -
# eine sichtbar durchscheinende Flaeche, hinter der nichts zu sehen ist.
# Also: erst das opake Geruest, dann die transparente Flaeche.
mol delrep 0 $espmol
set REP_MOL  0
set REP_SURF 1

_try color Name C gray
mol representation Licorice 0.150000 24.000000 24.000000
mol color Name
mol selection {all}
mol material Opaque
mol addrep $espmol

# Isoflaeche rho = ISO, eingefaerbt nach dem ZWEITEN Volumendatensatz. Das ist
# das VMD-Aequivalent zu ramp_new + set surface_color in PyMOL: VMD kennt kein
# Rampenobjekt, sondern faerbt direkt nach den Werten eines anderen Gitters in
# derselben Molekuel-ID.
#   Isosurface <iso> <volID> <show 0=Flaeche> <draw 0=solide> <step> <size>
mol representation Isosurface $ISO $VOL_RHO 0 0 1 1
mol color Volume $VOL_ESP
mol selection {all}
mol material Transparent
mol addrep $espmol

# @@RANGE_NEG@@ .. +@@RANGE@@ a.u. = @@KJ@@ kJ/(mol*e). RWB: rot negativ.
color scale method @@COLORSCALE@@
color scale midpoint 0.5
mol scaleminmax $espmol $REP_SURF @@RANGE_NEG@@ @@RANGE@@

# --- 3) Befehle -----------------------------------------------

# Zoom aus der Molekuelgroesse statt aus einer geratenen Zahl.
#
# Der sichtbare Weltausschnitt ist in der Hoehe rund "display height"/2 gross -
# empirisch an einem echten Tachyon-Render nachgemessen, VMD gibt die Groesse
# nirgends direkt aus. Damit das Molekuel FILL der Bildhoehe fuellt, muss also
# gelten:  scale = FILL * (height/2) / (2*r).
#
# r ist der Radius der Huellkugel plus rund 1.8 A fuer die Isoflaeche - ueber
# die Kugel, damit keine der drei Ansichten angeschnitten wird.
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
    # Notfallwert, falls diese VMD-Version die Hoehe nicht herausgibt: 6.0 ist
    # der Default. Ohne die Absicherung reisst ein leerer Rueckgabewert die
    # ganze Szene mit, weil esp_fit schon beim Start aufgerufen wird.
    if {[catch {set h [display get height]}] || ![string is double -strict $h]
        || $h <= 0} {
        set h 6.0
    }
    scale to [expr {$FILL * 0.25 * $h / $r}]
}

# Auf die ATOME zentrieren, nicht auf die Gitterbox: resetview passt die Ansicht
# an alle Reps an, und die Isoflaeche traegt die Ausdehnung des ganzen Gitters.
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

# Annahme: planares Molekuel in der xy-Ebene, Halogen bei y = 0, Ring bei y > 0
# (so legt die Turbomole-Optimierung die Halogenbenzole ab). Der Drehsinn bei
# sigma ist nicht beliebig - falsch herum zeigt das Bild die Gegenseite, und
# weil dort ebenfalls eine runde gefaerbte Flaeche sitzt, faellt das nicht auf.
proc esp_view {which} {
    global cbmol cbopts
    esp_center
    switch -- $which {
        pi      { }
        edge    { rotate y by 90 }
        sigma   { rotate x by -90
                  puts "Hinweis: in der Achsenansicht stapeln sich die"
                  puts "  transparenten Lagen - esp_opacity 1.0 oder"
                  puts "  esp_snapshot (Tachyon) gibt ein sauberes Bild." }
        default { puts "esp_view: pi | edge | sigma" ; return }
    }
    if {[info exists cbmol] && [lsearch [molinfo list] $cbmol] >= 0} {
        eval esp_colorbar $cbopts
    }
    puts "Ansicht: $which"
}

proc esp_iso {value} {
    global espmol VOL_RHO REP_SURF ISO
    set ISO $value
    mol modstyle $REP_SURF $espmol Isosurface $value $VOL_RHO 0 0 1 1
    puts "Isowert: rho = $value a.u."
}

proc esp_range {half} {
    global espmol REP_SURF RANGE cbmol
    set RANGE $half
    mol scaleminmax $espmol $REP_SURF [expr {-1.0 * $half}] $half
    if {[info exists cbmol] && [lsearch [molinfo list] $cbmol] >= 0} {
        esp_colorbar
    }
    puts "Farbskala: +/- $half a.u."
}

# Der Blick geht durch ZWEI Lagen der geschlossenen Flaeche; bei Deckkraft a
# bleibt (1-a)^2 uebrig. 0.5 ist der brauchbare Kompromiss, 1.0 zeigt nur die
# Flaeche, unter 0.3 wird die Farbe zu blass.
proc esp_opacity {value} {
    _try material change opacity Transparent $value
    puts "Deckkraft: $value"
}

# VMD hat keine Legende. Der Balken wird in eine eigene, leere Molekuel-ID
# gezeichnet und mit "mol fix" von der Maus abgekoppelt. Platziert wird er
# relativ zum Molekuel in Angstroem - in geratenen Bildschirmkoordinaten haengt
# das an Fenstergroesse und Zoom und landet lautlos ausserhalb des Bildes.
# Argumente: Abstand zur Huellkugel (A), Balkenhoehe (A), Schriftgroesse.
proc esp_colorbar {{gap 0.5} {barheight 0.7} {textsize 0.8}} {
    global espmol RANGE cbmol cbopts
    set cbopts [list $gap $barheight $textsize]

    set sel [atomselect $espmol all]
    set c  [measure center $sel]
    set mm [measure minmax $sel]
    $sel delete
    set ex [expr {[lindex [lindex $mm 1] 0] - [lindex [lindex $mm 0] 0]}]
    set ey [expr {[lindex [lindex $mm 1] 1] - [lindex [lindex $mm 0] 1]}]
    set ez [expr {[lindex [lindex $mm 1] 2] - [lindex [lindex $mm 0] 2]}]
    # Huellkugel statt Box, damit der Balken auch nach dem Drehen aussen bleibt.
    set r  [expr {0.5 * sqrt($ex*$ex + $ey*$ey + $ez*$ez) + 1.8}]
    set cx [lindex $c 0] ; set cy [lindex $c 1] ; set cz [lindex $c 2]
    set hw [expr {0.85 * $r}]
    set y0 [expr {$cy - $r - $gap}]
    set y1 [expr {$y0 + $barheight}]
    set x0 [expr {$cx - $hw}]

    if {[info exists cbmol] && [lsearch [molinfo list] $cbmol] >= 0} {
        graphics $cbmol delete all
    } else {
        mol new
        set cbmol [molinfo top]
        mol rename $cbmol colorbar
    }
    graphics $cbmol materials off   ;# Farben der Skala, nicht schattiert

    set c0 [colorinfo num]              ;# die Skala belegt die IDs oberhalb
    set c1 [expr {[colorinfo max] - 1}] ;# der benannten Farben
    set nseg 64
    set dx [expr {2.0 * $hw / $nseg}]
    for {set i 0} {$i < $nseg} {incr i} {
        graphics $cbmol color \
            [expr {int($c0 + double($i) / ($nseg - 1) * ($c1 - $c0))}]
        set xa [expr {$x0 + $i * $dx}]
        set xb [expr {$xa + $dx}]
        graphics $cbmol triangle "$xa $y0 $cz" "$xb $y0 $cz" "$xb $y1 $cz"
        graphics $cbmol triangle "$xa $y0 $cz" "$xb $y1 $cz" "$xa $y1 $cz"
    }

    graphics $cbmol color black
    set ty [expr {$y0 - 1.1 * $barheight}]
    graphics $cbmol text "$x0 $ty $cz" [format "%.3f" [expr {-1.0*$RANGE}]] \
        size $textsize thickness 2
    graphics $cbmol text "[expr {$cx - 0.3}] $ty $cz" "0" \
        size $textsize thickness 2
    graphics $cbmol text "[expr {$cx + 0.45*$hw}] $ty $cz" \
        [format "+%.3f" $RANGE] size $textsize thickness 2
    graphics $cbmol text "[expr {$cx - 1.2}] [expr {$y1 + 0.4*$barheight}] $cz" \
        "ESP / a.u." size $textsize thickness 2

    # Gleiche Zentrierung und Skalierung wie das Molekuel, aber ohne dessen
    # Rotation; mol fix koppelt zusaetzlich von der Maus ab.
    molinfo $cbmol set center_matrix [list [transoffset [vecscale -1.0 $c]]]
    molinfo $cbmol set rotate_matrix [list [transidentity]]
    molinfo $cbmol set global_matrix [molinfo $espmol get global_matrix]
    molinfo $cbmol set scale_matrix  [molinfo $espmol get scale_matrix]
    mol fix $cbmol
    mol top $espmol
}

proc esp_colorbar_off {} {
    global cbmol
    if {[info exists cbmol] && [lsearch [molinfo list] $cbmol] >= 0} {
        mol delete $cbmol
        unset cbmol
    }
}

# Strahlverfolgt. Schnelle Alternative: render snapshot $name.png
proc esp_snapshot {name} {
    render TachyonInternal $name.tga
    puts "geschrieben: $name.tga"
}

# --- 4) Anzeige -----------------------------------------------
_try display projection Orthographic
_try display depthcue off
_try display shadows off
_try display culling off
_try axes location off
_try color Display Background white
# Ohne das schneidet die nahe Clipping-Ebene beim Zoomen eine Scheibe heraus.
_try display nearclip set 0.010000
# Echte Transparenz statt Rasterpunkten ("screen door"); alte Treiber koennen
# kein GLSL - dann bleibt das Raster, die Tachyon-Bilder sind trotzdem sauber.
_try display rendermode GLSL
_try material change opacity  Transparent @@OPACITY@@
_try material change diffuse  Transparent 0.750000
# Glanzlicht klein halten - ein milchiger Schleier frisst den Durchblick auf.
_try material change specular Transparent 0.100000
_try material change shininess Transparent 0.300000

# --- 5) Start -------------------------------------------------
esp_view pi
if {[catch {esp_colorbar} err]} { puts "! Farbbalken: $err" }

puts ""
puts "esp.tcl geladen. rho = $ISO a.u., Farbskala +/- $RANGE a.u.@@STATS@@"
puts ""
