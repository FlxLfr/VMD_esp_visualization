# ==============================================================
# render_esp.tcl - die drei Standardansichten strahlverfolgt nach images/
#
# Aufruf aus dem Molekuelordner (dort, wo esp.tcl und die Cubes liegen):
#
#     vmd -e ../../scripts/render_esp.tcl
#
# Ueblicherweise ruft render_esp.py dieses Skript auf und macht danach die
# Farbskala und settings.txt. Direkt gestartet tut es nur den VMD-Teil.
#
# Vorher setzbar (in der Tk Console, dann "source"):
#   set ESP_RES    {1600 1280}   Fenstergroesse = Bildgroesse
#   set ESP_QUIT   0             VMD nach dem Rendern offen lassen
#   set ESP_AO     0             Umgebungsverdeckung aus
#   set ESP_OPAQUE 1             Isoflaeche opak rendern
#   set ESP_BG     black          Hintergrundfarbe
#   set ESP_VIEWS  {sigma}       nur einzelne Ansichten
#   set ESP_SNAPSHOT 1           Fenstermitschnitt statt Strahlverfolgung
#   set ESP_SCENE  esp_check.tcl andere Szenendatei (Selbsttest)
#   set ESP_OUTDIR images_check   anderer Zielordner (Selbsttest)
#   set ESP_PREFIX brombenzol     Dateipraefix (Standard: Ordnername)
# ==============================================================

if {![info exists ESP_RES]}    { set ESP_RES {1600 1280} }
if {![info exists ESP_QUIT]}   { set ESP_QUIT 1 }
if {![info exists ESP_AO]}     { set ESP_AO 1 }
if {![info exists ESP_OPAQUE]}  { set ESP_OPAQUE 0 }
if {![info exists ESP_BG]}     { set ESP_BG white }
# Angehaengt an den Dateinamen - bei mehreren Hintergrundfarben in einem Lauf
# unterscheidet er die Bildersaetze (<molekuel>_pi_black.tga).
if {![info exists ESP_SUFFIX]} { set ESP_SUFFIX "" }
if {![info exists ESP_VIEWS]}  { set ESP_VIEWS {pi edge sigma} }
if {![info exists ESP_SNAPSHOT]} { set ESP_SNAPSHOT 0 }
# Der Selbsttest schreibt seine Szene als esp_check.tcl, damit er die
# committete esp.tcl eines echten Laufs nicht ueberschreibt.
if {![info exists ESP_SCENE]}  { set ESP_SCENE esp.tcl }
# Zielordner und Praefix kommen von render_espVMD.py. Stehen sie nicht da,
# gelten dieselben Vorgaben wie beim Aufruf von Hand: images/ und der Name des
# Molekuelordners.
if {![info exists ESP_OUTDIR]} { set ESP_OUTDIR images }
if {![info exists ESP_PREFIX]} { set ESP_PREFIX [file tail [pwd]] }

if {![file exists $ESP_SCENE]} {
    puts "render_esp.tcl: $ESP_SCENE nicht gefunden."
    puts "  Aus dem Molekuelordner starten, oder erst xyzToCubeToVMDVis.py laufen lassen."
    if {$ESP_QUIT} { quit }
    return
}

source $ESP_SCENE

# Notausgang fuer die Achsenansicht: dort schneidet jeder Sehstrahl viele
# transparente Lagen der Isoflaeche, und mit Umgebungsverdeckung zusammen
# bringt das Tachyon zum Absturz. Opak gerendert faellt die Transparenz-
# rekursion weg. render_esp.py schaltet das nur zu, wenn es noetig war.
if {$ESP_OPAQUE} { esp_opacity 1.0 }

# --- Renderqualitaet ------------------------------------------
# Diese Einstellungen wirken im OpenGL-Fenster kaum und im Strahlverfolger
# stark. Vor allem die Umgebungsverdeckung: sie legt weiche Schatten in die
# Vertiefungen zwischen den CH-Ausbuchtungen und ist der groesste Einzelschritt
# von "Vorschau" zu "Publikationsbild".
# Die Umgebungsverdeckung wirft keine objektfoermigen Flecken, sondern dunkelt
# nur Vertiefungen ab - deshalb ist sie zuschaltbar.
#
# Schlagschatten dagegen bleiben fest aus, ohne Option. Sie werfen die
# Licorice-Staebchen als graue Kapseln auf die Isoflaeche: im pi-Bild sitzt
# hinter jedem Staebchen ein versetzter Doppelgaenger, in der Achsenansicht
# liegt ein Stab-Schatten mitten auf der Kugel. Beides sieht aus wie ein
# Artefakt der Daten und ist keines, und PyMOL rendert ebenfalls ohne - das
# haelt die Bilder vergleichbar. Wer sie doch sehen will, schaltet sie in der
# Tk Console zu: 'display shadows on', dann 'esp_snapshot <name>'.
_try display shadows off
_try color Display Background $ESP_BG
if {$ESP_AO} {
    _try display ambientocclusion on
    _try display aoambient 0.80
    _try display aodirect 0.40
}
_try display antialias on
_try display depthcue off
_try display resize [lindex $ESP_RES 0] [lindex $ESP_RES 1]
display update

# TachyonInternal rendert in Fenstergroesse. Wenn der Bildschirm kleiner ist
# als ESP_RES, klemmt Windows das Fenster ab und das Bild wird entsprechend
# kleiner - deshalb wird die tatsaechliche Groesse gemeldet.
puts "Fenster: [display get size]" ; flush stdout

# Der Fenstermitschnitt ist der Notausgang fuer die Achsenansicht: Tachyon
# steigt dort an der Zahl der durchquerten transparenten Lagen aus, OpenGL
# zeichnet sie dagegen ohne Rekursion. Das Bild ist etwas weniger fein, behaelt
# aber die Transparenz - besser als eine opake Ersatzdarstellung.
# Achtung: der Mitschnitt kopiert den Bildschirminhalt, das VMD-Fenster darf
# also nicht verdeckt sein.
#
# Sonst: auf der GPU liefert der OptiX-Pfad dieselbe Szene schneller. Nicht
# jeder Build hat ihn.
set RENDERER TachyonInternal
if {$ESP_SNAPSHOT} {
    set RENDERER snapshot
} elseif {[lsearch [render list] TachyonLOptiXInternal] >= 0} {
    set RENDERER TachyonLOptiXInternal
}
puts "Renderer: $RENDERER" ; flush stdout

# --- Rendern --------------------------------------------------
file mkdir $ESP_OUTDIR
set prefix $ESP_PREFIX
puts "Ziel: $ESP_OUTDIR" ; flush stdout

# Jede Ansicht einzeln absichern: bricht eine ab - Tachyon kann bei grossen
# Isoflaechen mit Umgebungsverdeckung aussteigen -, sollen die anderen trotzdem
# entstehen. Und nach jeder Zeile flush, damit im Log steht, wie weit VMD kam,
# falls der Prozess mittendrin stirbt.
foreach view $ESP_VIEWS {
    set out [file join $ESP_OUTDIR "${prefix}_${view}${ESP_SUFFIX}.tga"]
    puts "== $view: beginne ==" ; flush stdout
    if {[catch {esp_view $view ; display update} err]} {
        puts "! $view: Ansicht fehlgeschlagen: $err" ; flush stdout
        continue
    }
    if {[catch {render $RENDERER $out} err]} {
        puts "! $view: Rendern fehlgeschlagen: $err" ; flush stdout
    } else {
        puts "-> $out" ; flush stdout
    }
}

puts "render_esp.tcl fertig. PNG und Farbskala macht render_esp.py."
if {$ESP_QUIT} { quit }
