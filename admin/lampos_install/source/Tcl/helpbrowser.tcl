
set text ""
set labinfo ""
set yesno "no"
set altstart 1.0
set tags(none) ""
set maintopics 0
set topics ""
global srcpath; set srcpath .

proc gv_searchMatches {w pattern script} {
  scan [$w index end] %d numLines
  for {set i 1} {$i < $numLines} {incr i} {
    $w mark set last $i.0
    while {[regexp -nocase -indices $pattern [$w get last "last lineend"] indices]} {
      $w mark set first "last + [lindex $indices 0] chars"
      $w mark set last "last + 1 chars + [lindex $indices 1] chars"
      uplevel $script
      }
    }
  }
  
  proc gv_endofText {} {
    global yesno
    catch {destroy .eotw}
    toplevel .eotw -class Dialog
    wm withdraw .eotw
    wm title .eotw "Start at Top?"
    frame .eotw.main
    frame .eotw.bottom
    message .eotw.main.msg -text "You have reached the end of the text. Do you want to start at top?" -width 12c
    pack .eotw.main .eotw.bottom
    pack .eotw.main.msg -in .eotw.main
    button .eotw.bottom.yes -text "Yes" -command {
      set yesno yes
      destroy .eotw
      set start 1.0
      .hw.top.right.t1 yview 1.0
      catch {set index [.hw.top.right.t1 search -nocase -count count $searchvar $start]}
      .hw.top.right.t1 yview $index
      }
    button .eotw.bottom.no -text "No" -command {set yesno no; destroy .eotw}
    pack .eotw.bottom.yes .eotw.bottom.no -in .eotw.bottom -side left
    wm withdraw .eotw
    update idletasks
    set x [expr [winfo screenwidth .eotw]/2 -[winfo reqwidth .eotw]/2 - [winfo vrootx [winfo parent .eotw]]]
    set y [expr [winfo screenheight .eotw]/2 -[winfo reqheight .eotw]/2 - [winfo vrooty [winfo parent .eotw]]]
    wm geom .eotw +$x+$y
    wm deiconify .eotw
    grab set .eotw
  }

  proc gv_endofText2 {} {
    global yesno
    catch {destroy .eotw}
    toplevel .eotw -class Dialog
    wm withdraw .eotw
    wm title .eotw "Start at Top?"
    frame .eotw.main
    frame .eotw.bottom
    message .eotw.main.msg -text "You have reached the end of the text. Do you want to start at top?" -width 12c
    pack .eotw.main .eotw.bottom
    pack .eotw.main.msg -in .eotw.main
    button .eotw.bottom.yes -text "Yes" -command {
      set yesno yes
      destroy .eotw
      set start 1.0
      .licw.top.t1 yview 1.0
      catch {set index [.licw.top.t1 search -nocase -count count $searchvar $start]}
      .licw.top.t1 yview $index
      }
    button .eotw.bottom.no -text "No" -command {set yesno no; destroy .eotw}
    pack .eotw.bottom.yes .eotw.bottom.no -in .eotw.bottom -side left
    wm withdraw .eotw
    update idletasks
    set x [expr [winfo screenwidth .eotw]/2 -[winfo reqwidth .eotw]/2 - [winfo vrootx [winfo parent .eotw]]]
    set y [expr [winfo screenheight .eotw]/2 -[winfo reqheight .eotw]/2 - [winfo vrooty [winfo parent .eotw]]]
    wm geom .eotw +$x+$y
    wm deiconify .eotw
    grab set .eotw    
  }

  proc gv_notfound {} {
    global yesno
    catch {destroy .nfw}
    toplevel .nfw -class Dialog
    wm withdraw .nfw
    wm title .nfw "Start at Top?"
    frame .nfw.main
    frame .nfw.bottom
    message .nfw.main.msg -text "Search pattern not found!" -width 12c
    pack .nfw.main .nfw.bottom
    pack .nfw.main.msg -in .nfw.main
    button .nfw.bottom.yes -text "OK" -command {
      destroy .nfw}
    pack .nfw.bottom.yes -in .nfw.bottom -side left
    wm withdraw .nfw
    update idletasks
    set x [expr [winfo screenwidth .nfw]/2 -[winfo reqwidth .nfw]/2 - [winfo vrootx [winfo parent .nfw]]]
    set y [expr [winfo screenheight .nfw]/2 -[winfo reqheight .nfw]/2 - [winfo vrooty [winfo parent .nfw]]]
    wm geom .nfw +$x+$y
    wm deiconify .nfw
    grab set .nfw    
  }

###########################################################################
# Helpbrowser main routine:                                               #
###########################################################################

proc gv_helpbrowser {} {

global yesno labinfo text altstart maintopics topics srcpath

upvar tags tags

  catch {destroy .hw}
  toplevel .hw 
  wm title .hw Help
  wm transient .hw .

# frame .hw -height 15c -width 20c
# pack propagate .hw 0
# pack .hw
frame .hw.top
frame .hw.top.right
frame .hw.top.left
frame .hw.bottom

pack .hw.top .hw.bottom -in .hw -side top -expand yes
pack .hw.top.left .hw.top.right -in .hw.top -side left -expand yes

###########################################################################
# Define Listbox:                                                         #
###########################################################################

listbox .hw.top.left.content -relief sunken -borderwidth 6  -height 24 -width 27 -background white -foreground black -font -Adobe-Times-Medium-R-Normal--*-180-*
#scrollbar .hw.top.left.scroll -command ".hw.top.left.content yview"

pack .hw.top.left.content  -in .hw.top.left -fill y -side left -expand yes

###########################################################################
# Bind Listbox:                                                           #
###########################################################################

 bind .hw.top.left.content <Button-1> {
  after 100 {
  .hw.top.right.t1 yview $tags([selection get])
  }
  bind .hw.top.left.content <Double-Button-1> {
    after 100 {

      catch {set sel [selection get]}
      for {set maint 0} {$maint < $maintopics} {incr maint} {
        set s [lindex $topics $maint]
        if {$sel == [lindex $s 1]} {
          
          if {[lindex $s 0] == "+"} {                
            set s [lreplace $s 0 0 "-"]
            } else {
            set s [lreplace $s 0 0 "+"]
            }
          set topics [lreplace $topics $maint $maint $s]
          }
        }
        
    .hw.top.left.content delete 0 end
    for {set maint 0} {$maint < $maintopics} {incr maint} {
      set s [lindex $topics $maint]
      set labinfo [llength $s]
      if {[lindex $s 0]=="+"} {
        for {set submaint 1} {$submaint < [llength $s]} {incr submaint} {
          .hw.top.left.content insert end [lindex $s $submaint]}
        } else {
          .hw.top.left.content insert end [lindex $s 1]}
        }
      }  
    }
  }

###########################################################################
# Define Buttons/Labels/Entrys:                                           #
###########################################################################

label .hw.bottom.lab -text Searchpattern: -height 3
# label .hw.bottom.lab -textvariable labinfo
button .hw.bottom.ok -text OK -command {destroy .hw}
entry .hw.bottom.sentry -width 30 -relief sunken -bd 2 -textvariable searchvar
button .hw.bottom.sbutton -text Search -command {
  if {$searchvar != ""} {        
.hw.top.right.t1 tag delete found        
  set totalNoOfLines [lindex [.hw.top.right.t1 index end] 0]
  set indexStartOfPage [lindex [.hw.top.right.t1 yview] 0]
  set start [.hw.top.right.t1 index \
   "[expr int($totalNoOfLines * $indexStartOfPage)].0 + 2l"]
  .hw.top.right.t1 yview $start
 gv_searchMatches .hw.top.right.t1 $searchvar {
    .hw.top.right.t1 tag add found first last}
    .hw.top.right.t1 tag  configure found -background red

    set index [.hw.top.right.t1 search -nocase -count count $searchvar $start]
    set error [catch {.hw.top.right.t1 yview $index}]
    if {$error==1} {gv_notfound}
    set labinfo $altstart
 if {[.hw.top.right.t1 yview]==$altstart} {
       gv_endofText
}
    set altstart [.hw.top.right.t1 yview]
    } else {

.hw.top.right.t1 tag delete found
     }}

pack .hw.bottom.ok .hw.bottom.lab .hw.bottom.sentry .hw.bottom.sbutton -in .hw.bottom -side left 

###########################################################################
# Define Textfield:                                                       #
###########################################################################

  text .hw.top.right.t1 -height 28 -width 60 -yscrollcommand ".hw.top.right.scroll set" -wrap word -cursor arrow -font -Adobe-Times-Medium-R-Normal--*-180-* -background white -foreground black -borderwidth 3 -relief sunken
  scrollbar .hw.top.right.scroll -command ".hw.top.right.t1 yview"
  pack .hw.top.right.t1 .hw.top.right.scroll -in .hw.top.right -side left -fill y 

# Fill Texfield
  cd $srcpath
  set Datei [open "../../help.html" r]
  while {[gets $Datei Zeile] >= 0} {
  .hw.top.right.t1 insert end $Zeile\n}
  close $Datei

# Fix window size
wm resizable .hw 0 0


 set indexhl 0
 
###########################################################################
# Headlines:                                                              #
###########################################################################

 set maintopics 0
 gv_searchMatches .hw.top.right.t1 {(<h[0-9]>)([',0-9,".",a-z," ",A-Z,:,-,_,-]*)(</h[0-9]>)} {
   .hw.top.right.t1 tag add hline$indexhl first last
   set indices [.hw.top.right.t1 tag ranges hline$indexhl]   
   
   set headline [.hw.top.right.t1 get [lindex $indices 0] [lindex $indices 1]]
   set number $headline
   regsub {([',0-9,".",a-z," ",A-Z,:,-,_,-]*)(</h[0-9]>)} $number {} number
   regsub {<h} $number {} number
   regsub {>} $number {} number

   regsub {<h[0-9]>} $headline {} headline
   regsub {</h[0-9]>} $headline {} headline
   if {$number == 1} {
     .hw.top.right.t1 tag add mainhline first last
     }
   switch $number {
      1 {set headline "$headline"}
      2 {set headline "    $headline"}
      3 {set headline "        $headline"}
      4 {set headline "            $headline"}
      5 {set headline "                $headline"}
      6 {set headline "                    $headline"}
      7 {set headline "                        $headline"}
      }
      
   set tags($headline) [lindex $indices 0]

# Fill Listbox:
   if {$number == 1} {
     incr maintopics
     lappend maintopics$maintopics - $headline
     } else {
     lappend maintopics$maintopics $headline
     }
# .hw.top.left.content insert end $headline
   incr indexhl
   }
  for {set maint 1} {$maint <= $maintopics} {incr maint} {
    eval set s "\$maintopics$maint"
    lappend topics $s
    }
  for {set maint 0} {$maint < $maintopics} {incr maint} {
    set s [lindex $topics $maint]
    set labinfo [llength $s]
    if {[lindex $s 0]=="+"} {
    for {set submaint 1} {$submaint < [llength $s]} {incr submaint} {
      .hw.top.left.content insert end [lindex $s $submaint]}
      } else {
      .hw.top.left.content insert end [lindex $s 1]}
    }

  for {set index 0} {$index <= $indexhl} {incr index} {
    .hw.top.right.t1 tag configure hline$index \
    -font -adobe-times-Bold-R-Normal-*-*-180-*-*-*-*-*-*
    }

    .hw.top.right.t1 tag configure mainhline \
    -font -adobe-times-Bold-R-Normal-*-*-240-*-*-*-*-*-*

###########################################################################
# Deleting <>                                                             #
###########################################################################
    .hw.top.right.t1 tag configure italic \
    -font -adobe-times-medium-i-Normal-*-*-180-*-*-*-*-*-*

    .hw.top.right.t1 tag configure bold \
    -font -adobe-times-Bold-r-Normal-*-*-180-*-*-*-*-*-*

    .hw.top.right.t1 tag configure indent \
    -lmargin2 25


#bullets
  gv_searchMatches .hw.top.right.t1 {(^)(<li>)(.*)($)} {
    .hw.top.right.t1 tag add indent first last}

     gv_searchMatches .hw.top.right.t1 {<li>} {
    .hw.top.right.t1 delete first last
    .hw.top.right.t1 insert first "  -"    
    .hw.top.right.t1 tag add bold "first - 2 chars" first
     }

#bold
  gv_searchMatches .hw.top.right.t1 {(<b>)([',0-9,".",a-z," ",A-Z,:,-,_,-]*)(</b>)} {
    .hw.top.right.t1 tag add bold first last}
  gv_searchMatches .hw.top.right.t1 {<b>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {</b>} {
    .hw.top.right.t1 delete first last}

#italics
  gv_searchMatches .hw.top.right.t1 {(<i>)([',0-9,".",a-z," ",A-Z,:,-,_,-]*)(</i>)} {
    .hw.top.right.t1 tag add italic first last}
  gv_searchMatches .hw.top.right.t1 {<i>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {</i>} {
    .hw.top.right.t1 delete first last}

#definition list
  gv_searchMatches .hw.top.right.t1 {(<dt>)(.*)($)} {
    .hw.top.right.t1 tag add indent first last}
  gv_searchMatches .hw.top.right.t1 {<dd>} {
    .hw.top.right.t1 insert first "\n      "
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {<dt>} {
    .hw.top.right.t1 delete first last}

  gv_searchMatches .hw.top.right.t1 {<h[0-9]>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {</h[0-9]>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {<br>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {<p[0-9]>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {</p[0-9]>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {<p>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {</p>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {<pre>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {</pre>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {<ul>} {
    .hw.top.right.t1 delete first last}
  gv_searchMatches .hw.top.right.t1 {</ul>} {
    .hw.top.right.t1 delete first last}

.hw.top.right.t1 configure -state disabled
  grab set .hw
  wm withdraw .hw
  update idletasks
  set x [expr [winfo screenwidth .hw]/2 -[winfo reqwidth .hw]/2 - [winfo vrootx [winfo parent .hw]]]
  set y [expr [winfo screenheight .hw]/2 -[winfo reqheight .hw]/2 - [winfo vrooty [winfo parent .hw]]]
  wm geom .hw +$x+$y
  wm deiconify .hw
  tkwait window .hw
}

