#!/usr/bin/env k8

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

var c, min_dp = 5, min_dp_strand = 2, min_ab = 0.2, min_dp_bulk = 15, max_alt_dp_bulk = 0, max_lt = 0, is_hap_cell = false, min_mapq = 50;
var flt_win = 100, min_het_ab_bulk = 0.3, min_het_dp_bulk = 5, write_calls = false, fn_var = null;
while ((c = getopt(arguments, "ha:b:cv:")) != null) {
	if (c == 'h') is_hap_cell = true;
	else if (c == 'a') max_alt_dp_bulk = parseInt(getopt.arg);
	else if (c == 'b') min_dp_bulk = parseInt(getopt.arg);
	else if (c == 'c') write_calls = true;
	else if (c == 'v') fn_var = getopt.arg;
}

if (arguments.length - getopt.ind == 0) {
	print("Usage: plp-joint.js [options] <joint.vcf>");
	print("Options:");
	print("  -h        haploid mode");
	print("  -b INT    min bulk read depth [" + min_dp_bulk + "]");
	print("  -a INT    max bulk ALT read depth [" + max_alt_dp_bulk + "]");
	print("  -v FILE   common SNPs []");
	exit(1);
}

var file, buf = new Bytes();

function print_record(x, labels, calls, snv) {
	if (calls == null) {
		var hit = [];
		for (var j = 1; j < x.p.length; ++j) {
			if (x.p[j].flt) continue;
			if (x.p[j].alt_called) ++snv[j - 1];
			if (x.p[j].adf[1] + x.p[j].adr[1] == 0) continue;
			hit.push([labels[j-1], x.p[j].adf[1], x.p[j].adr[1]].join(":"));
		}
		print('MV', x.ctg, x.pos, x.ref, x.alt, x.p[0].adf[0]+x.p[0].adr[0], x.p[0].adf[1]+x.p[0].adr[1], hit.length, hit.join("\t"));
	} else {
		for (var j = 1; j < x.p.length; ++j) {
			var c;
			if (x.p[j].flt) c = '.';
			else if (x.p[j].adf[0] + x.p[j].adr[0] + x.p[j].adf[1] + x.p[j].adr[1] == 0) c = '.';
			else if (x.p[j].adf[0] + x.p[j].adr[0] > 0 && x.p[j].adf[1] + x.p[j].adr[1] > 0) c = '*';
			else if (x.p[j].adf[0] + x.p[j].adr[0] == 0) c = '+';
			else c = '-';
			calls[j-1].push(c);
		}
	}
}

var var_map = new Map();
if (fn_var != null) {
	warn('Reading sites to filter...');
	file = new File(fn_var);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0][0] == '#') continue;
		var_map.put(t[0] + ':' + t[1]);
	}
	file.close();
}

warn('Calling...');
file = arguments[getopt.ind] == '-'? new File() : new File(arguments[getopt.ind]);
var samples = [], last = [], calls = [], fn = [], snv = [], ado = [[], []], n_het_bulk = 0;
while (file.readline(buf) >= 0) {
	var m, t = buf.toString().split("\t");
	if (t[0] == '#CHROM') {
		for (var i = 10; i < t.length; ++i)
			samples.push(t[i].replace(/\.bam$/, ""));
		for (var i = 0; i < samples.length; ++i) {
			calls[i] = [];
			ado[0][i] = ado[1][i] = 0;
			fn[i] = snv[i] = 0;
		}
		if (!write_calls) print('SM', samples.join("\t"));
		continue;
	} else if (t[0][0] == '#') continue; // skip header
	t[1] = parseInt(t[1]);

	// skip bad sites: mapQ
	if ((m = /AMQ=([\d,]+)/.exec(t[7])) != null) {
		var s = m[1].split(","), flt = false;
		for (var j = 0; j < s.length; ++j)
			if (parseInt(s[j]) < min_mapq)
				flt = true;
		if (flt) continue;
	}

	// parse VCF (this part works with multiple ALT alleles)
	var p = [];
	for (var i = 9; i < t.length; ++i) {
		var s = t[i].split(":");
		var lt = s[3] == '.'? 0 : parseInt(s[3]);
		var adf = s[1].split(",");
		var adr = s[2].split(",");
		if (adf.length != adr.length) throw Error("Inconsistent VCF");
		var dp_ref = 0, dp_alt = 0;
		for (var j = 0; j < adf.length; ++j) {
			adf[j] = parseInt(adf[j]);
			adr[j] = parseInt(adr[j]);
			if (j == 0) dp_ref += adf[j] + adr[j];
			else dp_alt += adf[j] + adr[j];
		}
		if (i == 9) {
			p.push({ flt:false, adf:adf, adr:adr, lt:0 });
		} else {
			var flt = false;
			if (is_hap_cell && dp_alt > 0 && dp_ref > 0) flt = true;
			if (lt > max_lt) flt = true;
			if (flt) p.push({ flt:true });
			else p.push({ flt:false, adf:adf, adr:adr, lt:lt });
		}
	}

	// only consider biallelic substitutions (TODO: make it more general)
	if (t[3].length != 1 || t[4].length != 1) continue;

	// skip bad sites
	var dp_bulk = p[0].adf[0] + p[0].adf[1] + p[0].adr[0] + p[0].adr[1];
	if (dp_bulk < min_dp_bulk) continue;

	// count ADO and FN
	var is_het = false;
	if (p[0].adf[0] > 0 && p[0].adf[1] > 0 && p[0].adr[0] > 0 && p[0].adr[1] > 0 && p[0].adf[0] + p[0].adr[0] >= min_het_dp_bulk && p[0].adf[1] + p[0].adr[1] >= min_het_dp_bulk) {
		if ((p[0].adf[0] + p[0].adr[0]) / dp_bulk >= min_het_ab_bulk && (p[0].adf[1] + p[0].adr[1]) / dp_bulk >= min_het_ab_bulk)
			is_het = true;
	}
	if (is_het) {
		++n_het_bulk;
		for (var j = 1; j < p.length; ++j) {
			if (p[j].flt) {
				++ado[0][j - 1];
				++ado[1][j - 1];
				continue;
			}
			if (is_hap_cell) {
				if (p[j].adf[0] + p[j].adr[0] + p[j].adf[1] + p[j].adr[1] == 0)
					++ado[0][j - 1], ++ado[1][j - 1];
			} else {
				if (p[j].adf[0] + p[j].adr[0] == 0) ++ado[0][j - 1]; // ref allele dropped
				if (p[j].adf[1] + p[j].adr[1] == 0) ++ado[1][j - 1]; // alt allele dropped
			}
		}
	}

	// test if ALT is callable
	if (is_hap_cell) {
		for (var j = 1; j < p.length; ++j) {
			p[j].alt_called = false;
			if (p[j].flt) continue;
			if (p[j].adf[1] >= min_dp_strand && p[j].adr[1] >= min_dp_strand && p[j].adf[1] + p[j].adr[1] >= min_dp) {
				if (p[j].adf[0] + p[j].adr[0] == 0)
					p[j].alt_called = true;
			}
		}
	} else {
		for (var j = 1; j < p.length; ++j) {
			p[j].alt_called = false;
			if (p[j].flt) continue;
			if (p[j].adf[1] >= min_dp_strand && p[j].adr[1] >= min_dp_strand && p[j].adf[1] + p[j].adr[1] >= min_dp) {
				var ref = p[j].adf[0] + p[j].adr[0], alt = p[j].adf[1] + p[j].adr[1];
				if (alt >= (ref + alt) * min_ab)
					p[j].alt_called = true;
			}
		}
	}
	if (is_het)
		for (var j = 1; j < p.length; ++j)
			if (!p[j].alt_called) ++fn[j - 1];

	// test SNV
	if (p[0].adf[1] + p[0].adr[1] > max_alt_dp_bulk) continue;
	var alt_f = 0, alt_r = 0;
	for (var j = 1; j < p.length; ++j) {
		if (p[j].flt) continue;
		alt_f += p[j].adf[1];
		alt_r += p[j].adr[1];
	}
	if (alt_f < min_dp_strand || alt_r < min_dp_strand || alt_f + alt_r < min_dp) continue;

	// print
	while (last.length && (last[0].ctg != t[0] || last[0].pos + flt_win < t[1])) {
		var x = last.shift();
		if (!x.flt) print_record(x, samples, write_calls? calls : null, snv);
	}

	var flt_this = false;
	if (var_map && var_map.get(t[0] + ':' + t[1]) != null)
		flt_this = true;
	for (var j = 0; j < last.length; ++j) { // TODO: use a more sophisticated method
		flt_this = true;
		last[j].flt = true;
	}

	last.push({ flt:flt_this, ctg:t[0], pos:t[1], p:p, ref:t[3], alt:t[4] });
}
while (last.length) {
	var x = last.shift();
	if (!x.flt) print_record(x, samples, write_calls? calls : null, snv);
}
if (write_calls) {
	for (var i = 0; i < samples.length; ++i)
		print(samples[i], (ado[1][i] / n_het_bulk).toFixed(4), calls[i].join(""));
} else {
	var fnr = [], corr_snv = [];
	for (var i = 0; i < samples.length; ++i) {
		if (!is_hap_cell) fnr[i] = (fn[i] / n_het_bulk).toFixed(4);
		else fnr[i] = ((fn[i] - .5 * n_het_bulk) / (.5 * n_het_bulk)).toFixed(4);
		corr_snv[i] = (snv[i] / (1 - fnr[i])).toFixed(2);
	}
	print('FN', fnr.join("\t"));
	print('NM', snv.join("\t"));
	print('CN', corr_snv.join("\t"));
}

buf.destroy();
file.close();
