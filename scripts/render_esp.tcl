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
#   set ESP_SHADOWS 1            Schlagschatten an (siehe unten - meist nicht)
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
if {![info exists ESP_SHADOWS]} { set ESP_SHADOWS 0 }
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
# Ohne Fenster (-dispdev text) gibt es kein Display: "display resize" und der
# GLSL-Rendermodus laufen ins Leere, und die Bildgroesse kommt stattdessen von
# "-size" auf der Kommandozeile. Tachyon rendert trotzdem - er zeichnet aus dem
# Szenengraphen, nicht aus dem OpenGL-Puffer.
if {![info exists ESP_HEADLESS]} { set ESP_HEADLESS 0 }

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
# Schlagschatten sind hier bewusst AUS. Sie werfen die Licorice-Staebchen als
# graue Kapseln auf die Isoflaeche - im pi-Bild sitzt hinter jedem Staebchen ein
# versetzter Doppelgaenger, in der Achsenansicht liegt ein Stab-Schatten mitten
# auf der Kugel. Beides sieht aus wie ein Artefakt der Daten und ist keines.
# PyMOL rendert ebenfalls ohne Schatten, das haelt die Bilder vergleichbar.
#
# Die Umgebungsverdeckung bleibt an: sie ist der Teil, der Tiefe erzeugt, und
# sie wirft keine objektfoermigen Flecken, sondern dunkelt nur Vertiefungen ab.
_try display shadows [expr {$ESP_SHADOWS ? "on" : "off"}]
if {$ESP_AO} {
    _try display ambientocclusion on
    _try display aoambient 0.80
    _try display aodirect 0.40
}
_try display antialias on
_try display depthcue off
if {!$ESP_HEADLESS} {
    _try display resize [lindex $ESP_RES 0] [lindex $ESP_RES 1]
    display update
}

# TachyonInternal rendert in Fenstergroesse. Wenn der Bildschirm kleiner ist
# als ESP_RES, klemmt Windows das Fenster ab und das Bild wird entsprechend
# kleiner - deshalb wird die tatsaechliche Groesse gemeldet. Ohne Fenster
# entfaellt diese Deckelung: dort gilt "-size" und damit genau ESP_RES.
if {$ESP_HEADLESS} {
    puts "Fenster: [lindex $ESP_RES 0] [lindex $ESP_RES 1] (ohne Anzeige)"
} else {
    puts "Fenster: [display get size]"
}
flush stdout

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
if {$ESP_SNAPSHOT && $ESP_HEADLESS} {
    # Sollte nicht vorkommen - render_espVMD.py startet den Mitschnitt-
    # Durchgang immer mit Anzeige. Wer das Skript von Hand aufruft, soll
    # aber nicht ratlos vor einer leeren Datei stehen.
    puts "! Fenstermitschnitt braucht eine Anzeige - ohne Fenster nicht"
    puts "  moeglich. Es wird strahlverfolgt."
    set ESP_SNAPSHOT 0
}
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
    set out [file join $ESP_OUTDIR "${prefix}_${view}.tga"]
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
