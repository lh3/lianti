var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

var c, min_novo_dp = 5, min_bulk_dp = 15, min_bulk_var_dp = 5;
while ((c = getopt(arguments, "n:m:b:")) != null) {
	if (c == 'n') min_novo_dp = parseInt(getopt.arg);
	else if (c == 'm') min_bulk_var_dp = parseInt(getopt.arg);
	else if (c == 'b') min_bulk_dp = parseInt(getopt.arg);
}

if (getopt.ind == arguments.length) {
	print("Usage: k8 plp-diff.js [options] <input.vcf>");
	exit(1);
}

var file = arguments[getopt.ind] == "-"? new File() : new File(arguments[getopt.ind]);
var buf = new Bytes();

var n_bulk_var = 0, n_bulk_het = 0;
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	if (t[0].charAt(0) == '#') continue;
	if (t[3].length != 1 || t[4].length != 1) continue; // SNV ONLY!!!
	t[3] = t[3].toUpperCase();
	var u = t[9].split(/[:,]/);
	var v = t[10].split(/[:,]/);
	if (u.length < 5 || v.length < 5) continue; // something is wrong
	for (var i = 1; i <= 4; ++i)
		u[i] = parseInt(u[i]), v[i] = parseInt(v[i]);
	if (t[3] > t[4]) {
		if (t[3] == 'C') t[3] = 'G';
		else if (t[3] == 'G') t[3] = 'C';
		else if (t[3] == 'T') t[3] = 'A';
		if (t[4] == 'A') t[4] = 'T';
		else if (t[4] == 'C') t[4] = 'G';
		else if (t[4] == 'G') t[4] = 'C';
	}
	if (u[1] + u[2] + u[3] + u[4] < min_bulk_dp) continue; // bulk does not have enough coverage
	++n_bulk_var;
	var tag = [];
	if (u[2] + u[4] == 0) { // non variant in bulk
		if (v[2] + v[4] >= min_novo_dp) tag.push('FP');
	} else if (u[2] + u[4] >= min_bulk_var_dp) {
		if (v[2] + v[4] == 0) tag.push('FN');
		if (u[1] + u[3] >= min_bulk_var_dp) {
			++n_bulk_het;
			if (v[1] + v[3] == 0) tag.push('ADO_R');
			if (v[2] + v[4] == 0) tag.push('ADO_A');
		}
	}
	if (tag.length) print(t[0], t[1], t[3], t[4], t[3]+t[4], tag.join(","), u[1] + u[3], u[2] + u[4], v[1] + v[3], v[2] + v[4], v.length > 5? v[5] : 0, t[7]);
}
warn(n_bulk_var, n_bulk_het);

buf.destroy();
file.close();
