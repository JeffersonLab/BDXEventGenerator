BEGIN {
  print "<html>"
  
  print "<head>"
  print "  <meta HTTP-EQUIV=\"Content-Type\" CONTENT=\"text/html; charset=iso-8859-1\">"
  print "  <meta NAME=\"keywords\" CONTENT=\"root, tree, ntuple, format, description\">"
  print "  <title>root tree description</title>"
  print "</head>"

  print "<body>"
  
  print "<H1>root tree description</H1>"
  print "<hr>"
  print "<H2>Branches</H2>"
  print "<hr>"

  print "<table style=\"border: 1px dotted;\" align=\"center\" border=\"0\" cellpadding=\"7\" cellspacing=\"3\" widt=\"95%\">"
  print "<tr><td><b>Branch</b></td>"
  print "<td><b>Definition</b></td>"
  print "<td><b>Class</b></td></tr>"

  even = 1;
  previous_line = ""
}

function print_line(previous_line, current_line, even) {
  if(even) print "<tr bgcolor=\"#eeeeee\">"
 	else print "<tr bgcolor=\"#ffffff\">"
  split(current_line, a, "\"");
  print "  <td>"a[2]"</td>"
  split(previous_line, a, "// ");
  print "  <td>"a[2]"</td>"
  split(current_line, a, ",");
  split(a[2], b, "::");
  gsub(" ", "", b[1]);
 
  if(b[1] == "TRootGen" || b[1] == "TRootSimParticle") b[1] = "TRootGenParticle";
 
  print "  <td><a href=\"#"b[1]"\">"b[1]"</a></td>";
  print "</tr>";
}

/NewBranch/ {
  print_line(previous_line, $0, even);
  even = !even;
}

{
  previous_line = $0
}

END {
  print "</table>"
}
