#!/usr/bin/env k8

/************
 * getopt() *
 ************/

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

/*************************************
 * Parameters & command-line parsing *
 *************************************/

var c, min_mapq = 50, flt_win = 100, n_bulk = 1, is_hap_cell = false;
var min_dp_alt_cell = 5, min_dp_alt_strand_cell = 2, min_ab_cell = 0.2, max_lt_cell = 0;
var min_dp_bulk = 15, max_alt_dp_bulk = 0, min_het_ab_bulk = 0.3, min_het_dp_bulk = 5;
var fn_var = null, fn_hap = null, fn_excl = null;

while ((c = getopt(arguments, "h:A:b:v:D:e:Hl:a:s:w:")) != null) {
	if (c == 'b') n_bulk = parseInt(getopt.arg);
	else if (c == 'H') is_hap_cell = true;
	else if (c == 'h') fn_hap = getopt.arg;
	else if (c == 'e') fn_excl = getopt.arg;
	else if (c == 'v') fn_var = getopt.arg;
	else if (c == 'a') min_dp_alt_cell = parseInt(getopt.arg);
	else if (c == 's') min_dp_alt_strand_cell = parseInt(getopt.arg);
	else if (c == 'w') flt_win = parseInt(getopt.arg);
	else if (c == 'l') max_lt_cell = parseInt(getopt.arg);
	else if (c == 'D') min_dp_bulk = parseInt(getopt.arg);
	else if (c == 'A') max_alt_dp_bulk = parseInt(getopt.arg);
}

if (arguments.length - getopt.ind == 0) {
	print("Usage: plp-joint.js [options] <joint.vcf>");
	print("Options:");
	print("  General:");
	print("    -b INT    number of bulk samples [1]");
	print("    -h FILE   samples in FILE are haploid []");
	print("    -H        mark all single-cell samples as haploid");
	print("    -e FILE   exclude samples contained in FILE []");
	print("    -v FILE   exclude positions in VCF FILE []");
	print("  Cell:");
	print("    -a INT    min ALT read depth to call an SNV [" + min_dp_alt_cell + "]");
	print("    -s INT    min ALT read depth per strand [" + min_dp_alt_strand_cell + "]");
	print("    -l INT    max LIANTI conflicting reads [" + max_lt_cell + "]");
	print("    -w INT    size of window to filter clustered SNVs [" + flt_win + "]");
	print("  Bulk:");
	print("    -D INT    min bulk read depth [" + min_dp_bulk + "]");
	print("    -A INT    max bulk ALT read depth [" + max_alt_dp_bulk + "]");
	exit(1);
}

/***********************
 * Auxiliary functions *
 ***********************/

function read_list(fn)
{
	if (fn == null || fn == "") return {};
	var buf = new Bytes();
	var file = fn == '-'? new File() : new File(fn);
	var h = {};
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		h[t[0]] = 1;
	}
	file.close();
	buf.destroy();
	return h;
}

function aggregate_calls(x, cell_meta)
{
	var bulk_ad = [0, 0], cell_hit = [];
	for (var i = 0; i < x.bulk.length; ++i)
		bulk_ad[0] += x.bulk[i].ad[0], bulk_ad[1] += x.bulk[i].ad[1];
	for (var i = 0; i < x.cell.length; ++i) {
		var c = x.cell[i], b;
		if (c.flt) b = '.';
		else if (c.dp == 0) b = '.';
		else if (c.ad[0] > 0 && c.ad[1] > 0) b = '*';
		else if (c.ad[1] > 0) b = '+';
		else b = '-';
		cell_meta[i].calls.push(b);
		if (!c.flt) {
			if (c.alt) ++cell_meta[i].snv;
			if (c.ad[1] > 0)
				cell_hit.push([cell_meta[i].name, c.adf[1], c.adr[1]].join(":"));
		}
	}
	print('NV', x.ctg, x.pos, x.ref, x.alt, bulk_ad.join("\t"), cell_hit.length, cell_hit.join("\t"));
}

/********
 * Main *
 ********/

var file, buf = new Bytes();

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

var sample_excl = read_list(fn_excl);
var sample_hap = read_list(fn_hap);
var col2cell = [];
var cell_meta = [];

warn('Calling...');
file = arguments[getopt.ind] == '-'? new File() : new File(arguments[getopt.ind]);
var last = [], n_het_bulk = 0;
while (file.readline(buf) >= 0) {
	var m, t = buf.toString().split("\t");
	if (t[0] == '#CHROM') { // parse the sample line
		var sample_name = [];
		for (var i = 9 + n_bulk; i < t.length; ++i) {
			var s1 = t[i], s2 = s1.replace(/\.bam$/, "");
			if (sample_excl[s1] || sample_excl[s2]) continue;
			var pl = sample_hap[s1] || sample_hap[s2]? 1 : 2;
			cell_meta.push({ name:s2, ploidy:pl, col:i, ado:[0,0], fn:0, snv:0, calls:[] });
			sample_name.push(s2); // for printing only
		}
		for (var i = 0; i < cell_meta.length; ++i)
			col2cell[cell_meta[i].col] = i;
		print('SM', sample_name.join("\t"));
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
	var cell = [], bulk = [];
	for (var i = 9; i < t.length; ++i) {
		var cell_id = col2cell[i];
		if (i >= 9 + n_bulk && cell_id == null) continue; // exclude this sample
		var s = t[i].split(":");
		var lt = s[3] == '.'? 0 : parseInt(s[3]);
		var adf = s[1].split(",");
		var adr = s[2].split(",");
		var ad = [], dp = 0;
		if (adf.length != adr.length) throw Error("Inconsistent VCF");
		var dp_ref = 0, dp_alt = 0;
		for (var j = 0; j < adf.length; ++j) {
			adf[j] = parseInt(adf[j]);
			adr[j] = parseInt(adr[j]);
			if (j == 0) dp_ref += adf[j] + adr[j];
			else dp_alt += adf[j] + adr[j];
			ad[j] = adf[j] + adr[j];
			dp += ad[j];
		}
		if (i < 9 + n_bulk) {
			bulk.push({ dp:dp, ad:ad, adf:adf, adr:adr });
		} else {
			var flt = false;
			if (cell_meta[cell_id].ploidy == 1 && dp_alt > 0 && dp_ref > 0) flt = true; // two alleles in a haploid cell
			if (lt > max_lt_cell) flt = true;
			cell.push({ flt:flt, dp:dp, ad:ad, adf:adf, adr:adr, lt:lt });
		}
	}

	// only consider biallelic substitutions (TODO: make it more general)
	if (t[3].length != 1 || t[4].length != 1) continue;

	// test het in the bulk(s)
	var bulk_dp_low = false, all_het = true;
	for (var i = 0; i < bulk.length; ++i) {
		var b = bulk[i];
		b.het = false;
		if (b.adf[0] > 0 && b.adf[1] > 0 && b.adr[0] > 0 && b.adr[1] > 0 && b.ad[0] >= min_het_dp_bulk && b.ad[1] >= min_het_dp_bulk) {
			if (b.ad[0] >= b.dp * min_het_ab_bulk && b.ad[1] >= b.dp * min_het_ab_bulk)
				b.het = true;
		}
		if (!b.het) all_het = false;
		if (b.dp < min_dp_bulk)
			bulk_dp_low = true;
	}
	if (bulk_dp_low) continue; // dp of one bulk is too low. Skip the rest

	// count ADO
	if (all_het) {
		++n_het_bulk;
		for (var j = 0; j < cell.length; ++j) {
			var c = cell[j];
			if (c.flt) { // if filtered, both alleles are dropped
				++cell_meta[j].ado[0];
				++cell_meta[j].ado[1];
				continue;
			}
			if (c.ad[0] == 0) ++cell_meta[j].ado[0]; // ref allele dropped
			if (c.ad[1] == 0) ++cell_meta[j].ado[1]; // alt allele dropped
		}
	}

	// test if ALT is callable and count FN
	for (var i = 0; i < cell.length; ++i) {
		var c = cell[i];
		// If a cell is haploid and it has ref alleles, c.flt will be true. The conditions below work with haploid cells.
		c.alt = (!c.flt && c.ad[1] >= min_dp_alt_cell && c.adf[1] >= min_dp_alt_strand_cell && c.adr[1] >= min_dp_alt_strand_cell && c.ad[1] >= c.dp * min_ab_cell);
		if (all_het && !c.alt) ++cell_meta[i].fn;
	}

	// test SNV
	var n_bulk_ref = 0;
	for (var i = 0; i < bulk.length; ++i)
		if (bulk[i].ad[1] <= max_alt_dp_bulk)
			++n_bulk_ref;
	if (n_bulk_ref == 0) continue;

	var cell_alt_f = 0, cell_alt_r = 0;
	for (var i = 0; i < cell.length; ++i) {
		if (cell[i].flt) continue;
		cell_alt_f += cell[i].adf[1];
		cell_alt_r += cell[i].adr[1];
	}
	if (cell_alt_f < min_dp_alt_strand_cell || cell_alt_r < min_dp_alt_strand_cell || cell_alt_f + cell_alt_r < min_dp_alt_cell) // too few ALT reads in cell(s)
		continue;

	// filter by window & print
	while (last.length && (last[0].ctg != t[0] || last[0].pos + flt_win < t[1])) {
		var x = last.shift();
		if (!x.flt) aggregate_calls(x, cell_meta);
	}

	var flt_this = false;
	if (var_map && var_map.get(t[0] + ':' + t[1]) != null)
		flt_this = true;
	for (var j = 0; j < last.length; ++j) { // TODO: use a more sophisticated method
		flt_this = true;
		last[j].flt = true;
	}

	last.push({ flt:flt_this, ctg:t[0], pos:t[1], bulk:bulk, cell:cell, ref:t[3], alt:t[4] });
}
while (last.length) {
	var x = last.shift();
	if (!x.flt) aggregate_calls(x, cell_meta);
}

/***************************
 * Output final statistics *
 ***************************/

var snv = [], fnr = [], corr_snv = [];
for (var i = 0; i < cell_meta.length; ++i) {
	var c = cell_meta[i];
	snv[i] = c.snv;
	if (c.ploidy == 2) fnr[i] = (c.fn / n_het_bulk).toFixed(4);
	else fnr[i] = ((c.fn - .5 * n_het_bulk) / (.5 * n_het_bulk)).toFixed(4);
	corr_snv[i] = (c.snv / (1.0 - fnr[i])).toFixed(2);
}
print('NN', snv.join("\t"));
print('FN', fnr.join("\t"));
print('CN', corr_snv.join("\t"));

// output "multi-alignment"
for (var i = 0; i < cell_meta.length; ++i)
	print('CA', cell_meta[i].name, (cell_meta[i].ado[1] / n_het_bulk).toFixed(4), cell_meta[i].calls.join(""));

/********
 * Free *
 ********/

if (var_map != null) var_map.destroy();
buf.destroy();
file.close();