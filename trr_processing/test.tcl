proc test {seltext} {

set nf [molinfo top get numframes]
set sel [atomselect top "$seltext"]

for {set i 1} {$i < $nf} {incr i} {
	$sel frame $i

	set cor1 [$sel get {x}]
	set cor2 [$sel get {y}]
	set cor3 [$sel get {z}]	

	puts "frame $i distance $cor1 $cor2 $cor3."

}

}
