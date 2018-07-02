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

var c, min_snv_dp = 5, min_snv_dp_ds = 1, min_snv_ab = .2, min_snv_frag_conflict = 0, min_bulk_dp = 15, min_bulk_var_dp = 5, min_bulk_het_ab = .3, min_mapq = 40, min_snv_dist = 100, max_bulk_cnt = 0;
var hap = false, max_hap_err = 1, output_TP = false, force_sgl = false, pair_mode = false, cnt_gap = false;
while ((c = getopt(arguments, "n:m:b:q:a:A:d:he:P1pgs:f:c:")) != null) {
	if (c == 'n') min_snv_dp = parseInt(getopt.arg);
	else if (c == 'c') max_bulk_cnt = parseInt(getopt.arg);
	else if (c == 's') min_snv_dp_ds = parseInt(getopt.arg);
	else if (c == 'm') min_bulk_var_dp = parseInt(getopt.arg);
	else if (c == 'b') min_bulk_dp = parseInt(getopt.arg);
	else if (c == 'q') min_mapq = parseInt(getopt.arg);
	else if (c == 'a') min_snv_ab = parseFloat(getopt.arg);
	else if (c == 'A') min_bulk_het_ab = parseFloat(getopt.arg);
	else if (c == 'd') min_snv_dist = parseInt(getopt.arg);
	else if (c == 'f') min_snv_frag_conflict = parseInt(getopt.arg);
	else if (c == 'h') hap = true;
	else if (c == 'P') output_TP = true;
	else if (c == 'e') max_hap_err = parseInt(getopt.arg);
	else if (c == '1') force_sgl = true;
	else if (c == 'p') pair_mode = true;
	else if (c == 'g') cnt_gap = true;
}

if (getopt.ind == arguments.length) {
	print("Usage: k8 plp-diff.js [options] <input.vcf>");
	print("Options:");
	print("  -q INT     min RMS mapping quality ["+min_mapq+"]");
	print("  -b INT     min bulk read depth to call a het ["+min_bulk_dp+"]");
	print("  -m INT     min bulk allele depth to call a het ["+min_bulk_var_dp+"]");
	print("  -A FLOAT   min bulk allele balance to call a het ["+min_bulk_het_ab+"]");
	print("  -n INT     min single-cell ALT read depth to call a SNV ["+min_snv_dp+"]");
	print("  -s INT     min single-cell ALT read depth on each strand ["+min_snv_dp_ds+"]");
	print("  -a FLOAT   min single-cell ALT allele balance to call a SNV ["+min_snv_ab+"]");
	print("  -f INT     min single-cell fragment conflicts ["+min_snv_frag_conflict+"]");
	print("  -d INT     drop SNVs within INT-bp between each other ["+min_snv_dist+"]");
	print("  -h         haploid mode");
	print("  -e INT     ignore a bulk variant if #ref_alleles > INT ["+max_hap_err+"]");
	print("  -1         only look at the first two samples");
	print("  -p         double-cell mode");
	print("  -g         count gaps");
	exit(1);
}

var file = arguments[getopt.ind] == "-"? new File() : new File(arguments[getopt.ind]);
var buf = new Bytes();

var n_bulk_het = 0, n_ado_ref = 0, n_ado_alt = 0, n_ado_both = 0, n_het_fn = 0, n_het_fn2 = 0, n_snv = 0, n_snv_nonCT = 0, n_ins = 0, n_del = 0;
var last = [];
while (file.readline(buf) >= 0) {
	var m, is_indel, t = buf.toString().split("\t");
	if (t[0].charAt(0) == '#') continue; // skip VCF header
	if (t.length < 11) {
		warn("WARNING: incomplete line: '" + buf.toString() + "'");
		continue;
	}
	is_indel = (t[3].length == 1 && t[4].length == 1)? false : true;
	if (!cnt_gap && is_indel) continue;
	if (force_sgl) t.length = 11;
	t[1] = parseInt(t[1]);
	t[3] = t[3].toUpperCase();
	var u = t[9].split(/[:,]/);
	var v = t[10].split(/[:,]/);
	var w = t.length >= 12? t[11].split(/[:,]/) : null;
	if (u.length < 5 || v.length < 5) continue; // something is wrong
	for (var i = 1; i < u.length; ++i) { // convert to integers
		u[i] = u[i] == '.'? 0 : parseInt(u[i]);
		v[i] = v[i] == '.'? 0 : parseInt(v[i]);
		if (w) w[i] = w[i] == '.'? 0 : parseInt(w[i]);
	}
	var ref = t[3], alt = t[4];
	if (t[3] > t[4] && !is_indel) { // determine mutation type
		if (t[3] == 'C') t[3] = 'G';
		else if (t[3] == 'G') t[3] = 'C';
		else if (t[3] == 'T') t[3] = 'A';
		if (t[4] == 'A') t[4] = 'T';
		else if (t[4] == 'C') t[4] = 'G';
		else if (t[4] == 'G') t[4] = 'C';
	}
	var bulk_dp = u[1] + u[2] + u[3] + u[4];
	if (bulk_dp < min_bulk_dp) continue; // bulk does not have enough coverage
	if ((m = /\bAMQ=([\.\d,]+)/.exec(t[7])) != null) { // actually for bialliac SNVs, the block can be simpler; but let's be more general
		var mq = 256, s = m[1].split(",");
		for (var i = 0; i < s.length; ++i) {
			if (s[i] == '.') continue;
			var x = parseInt(s[i]);
			mq = mq < x? mq : x;
		}
		if (mq < min_mapq) continue;
	}
	var is_snv_called;
	if (t.length <= 11) { // only two samples: bulk and a single cell
		if (!hap) is_snv_called = v[2] + v[4] >= min_snv_dp && (v[2] + v[4]) / (v[1] + v[2] + v[3] + v[4]) >= min_snv_ab? true : false;
		else is_snv_called = v[2] + v[4] >= min_snv_dp && v[1] + v[3] == 0? true : false;
		if (v.length > 5 && v[5] > min_snv_frag_conflict) is_snv_called = false;
		if (v[2] < min_snv_dp_ds || v[4] < min_snv_dp_ds) is_snv_called = false;
	} else { // three samples: bulk and a pair of single cells; FIXME: min_snv_dp_ds is not used in this block
		if (pair_mode && !hap) {
			is_snv_called = true;
			if (v[2] + v[4] + w[2] + w[4] < min_snv_dp) is_snv_called = false;
			else if ((v[2] + v[4] + w[2] + w[4]) / (v[1] + v[2] + v[3] + v[4] + w[1] + w[2] + w[3] + w[4]) < min_snv_ab) is_snv_called = false;
			else if (v[2] + v[4] < 3 || w[2] + w[4] < 3) is_snv_called = false;
		} else {
			var called1, called2;
			if (!hap) {
				called1 = v[2] + v[4] >= min_snv_dp && (v[2] + v[4]) / (v[1] + v[2] + v[3] + v[4]) >= min_snv_ab? true : false;
				called2 = w[2] + w[4] >= min_snv_dp && (w[2] + w[4]) / (w[1] + w[2] + w[3] + w[4]) >= min_snv_ab? true : false;
			} else {
				called1 = v[2] + v[4] >= min_snv_dp && v[1] + v[3] == 0? true : false;
				called2 = w[2] + w[4] >= min_snv_dp && w[1] + w[3] == 0? true : false;
			}
			is_snv_called = called1 && called2? true : false;
		}
		if (v.length > 5 && v[5] > 0) is_snv_called = false;
		if (w.length > 5 && w[5] > 0) is_snv_called = false;
	}
	if (!hap) {
		if (u[1] > 0 && u[2] > 0 && u[3] > 0 && u[4] > 0 && u[1] + u[3] >= min_bulk_var_dp && u[2] + u[4] >= min_bulk_var_dp
			&& (u[1] + u[3]) / bulk_dp >= min_bulk_het_ab && (u[2] + u[4]) / bulk_dp >= min_bulk_het_ab) // a bulk het
		{
			++n_bulk_het;
			// count ADO
			if (v[1] + v[2] + v[3] + v[4] == 0) ++n_ado_ref, ++n_ado_alt, ++n_ado_both;
			else if (v[1] + v[3] == 0) ++n_ado_ref;
			else if (v[2] + v[4] == 0) ++n_ado_alt;
			// count FN
			if (!is_snv_called) ++n_het_fn;
			if (is_snv_called && output_TP) print("TP", t[0], t[1], ref, alt, t[3]+t[4], u[1]+u[3], u[2]+u[4]);
		}
	} else {
		if (u[2] > 0 && u[4] > 0 && u[2] + u[4] >= min_bulk_var_dp && u[1] + u[3] <= max_hap_err) {
			++n_bulk_het; // ok, this is really a hom
			if (v[2] + v[4] == 0) ++n_ado_alt;
			if (!is_snv_called) ++n_het_fn;
		}
	}
	if (u[2] + u[4] <= max_bulk_cnt && is_snv_called) { // a potential SNV
		var s, type, flt_this = false;
		type = t[3].length > t[4].length? 'DEL' : t[3].length < t[4].length? 'INS' : t[3]+t[4];
		if (t.length <= 11) s = [t[0], t[1], ref, alt, type, v[1] + v[3], v[2] + v[4], v.length > 5? v[5] : 0, t[7]];
		else s = [t[0], t[1], ref, alt, type, v[1] + v[3] + w[1] + w[3], v[2] + v[4] + w[2] + w[4], v.length > 5? v[5] + w[5] : 0, t[7]];
		for (var i = 0; i < last.length; ++i) {
			if (last[i][0] != t[0] || t[1] - last[i][1] > min_snv_dist) {
				if (!last[i][2]) {
					print('NV', last[i][3].join("\t"));
					if (type == 'INS') ++n_ins;
					else if (type == 'DEL') ++n_del;
					else ++n_snv;
					if (last[i][3][4] != "CT" && type != 'INS' && type != 'DEL')
						++n_snv_nonCT;
				}
				last.shift();
				--i;
			} else last[i][2] = flt_this = true; // filtered
		}
		last.push([t[0], t[1], flt_this, s]);
	}
}
for (var i = 0; i < last.length; ++i)
	if (!last[i][2]) {
		print('NV', last[i][3].join("\t"));
		++n_snv;
		if (last[i][3][4] != "CT") ++n_snv_nonCT;
	}

print("NE", n_bulk_het);
print("VN", n_het_fn, (n_het_fn / n_bulk_het).toFixed(4));
if (!hap) print("RO", n_ado_ref, (n_ado_ref / n_bulk_het).toFixed(4));
print("AO", n_ado_alt, (n_ado_alt / n_bulk_het).toFixed(4));
if (!hap) print("BO", n_ado_both, (n_ado_both / n_bulk_het).toFixed(4));
print("NN", n_snv, n_snv_nonCT, n_del, n_ins);

buf.destroy();
file.close();
