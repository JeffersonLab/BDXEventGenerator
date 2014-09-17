BEGIN {  
  print "<hr>"
  print "<H2>Classes</H2>"
  print "<hr>"

  print "<table style=\"border: 1px dotted;\" align=\"center\" border=\"0\" cellpadding=\"7\" cellspacing=\"3\" widt=\"95%\">"
  print "<tr><td><b>Parameter</b></td>"
  print "<td><b>Definition</b></td>"
  print "<td><b>How it was calculated</b></td></tr>"

}

function print_line(name, comment, even, end) {
  if(name != ""){
   	if(even) print "<tr bgcolor=\"#eeeeee\">"
   	else print "<tr bgcolor=\"#ffffff\">"
    print "  <td>"name"</td>"
    split(comment, a, "|");
    print "  <td>"a[1]"</td>"
    print "  <td>"a[2]"</td>"
    print "</tr>"
  }
}

/^ *class TRoot/{
  print_line(name, comment, even, 1);
  even = 1;
  name = "";
  comment = "";
  split($2, a, ":");
  print "<tr bgcolor=\"#ffffff\"><td colspan=3><hr><a name=\""a[1]"\"><H3>"a[1]"</H3><hr></td></tr>"
}

/: public TRoot/{
  name = sprintf("<a href=\"#%s\">%s</a>", $4, $4);
  split($2, a, ":");
  comment = sprintf("%s inherits all %s parameters", a[1], $4);
}

/^ *[A-Za-z_]* [A-Za-z].*; \/\/ / {
  print_line(name, comment, even, 0);
  split($2, a, ";");
  name = a[1];
  split($0, a, "// ");
  comment = a[2];
  even = !even;
}

/^ +\/\/ /{split($0, a, "// "); comment = comment" "a[2]}
END {
  print_line(name, comment, even, 1);
  print "</table>"
  print "</body></html>"
}

