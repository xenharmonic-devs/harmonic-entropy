// import dsp.js

var plot;
var r = new Array();
var oldN, oldseries;
var plotted = false;

//Copyright 2009 Nicholas C. Zakas. All rights reserved.
//MIT Licensed

function conv(olda, oldb) {
	if(olda.length != oldb.length) {
		throw "conv(...): input arrays must be the same length";
		return null;
	}
		
	var a = olda.concat();
	var b = oldb.concat();
	
	var len = a.length;
	var minlen = 1;
	while (minlen	< len)	//2 is to avoid circular convolution distortion
		minlen *= 2;
	
	var padding = minlen - len;
	for(var i=len; i<minlen; i++) {
		a[i] = 0;
		b[i] = 0;
	}

	var f_a = new FFT(minlen,1);
	var f_b = new FFT(minlen,1);
	f_a.forward(a);
	f_b.forward(b);
	
	//(a+bi)(c+di)
	//(ac - bd) + (ad + bc)i
	
	//Do the multiplication out
	var f_out = new FFT(minlen,1);
	for(var i=minlen-1; i>=0; i--) {
		f_out.real[i] = f_a.real[i]*f_b.real[i] - f_a.imag[i]*f_b.imag[i];
		f_out.imag[i] = f_a.real[i]*f_b.imag[i] + f_a.imag[i]*f_b.real[i];
	}

	var out = f_out.inverse().slice(0,len);
	return out;
}

function log2(x) {
	return Math.log(x)/Math.log(2);
}

function gcd(x, y) {
	while (y != 0) {
		var z = x % y;
		x = y;
		y = z;
	}
	return x;
}

function harmonicEntropy(HEinfo) {
	//globals
	var a = parseFloat(HEinfo.a);
	var s = parseFloat(HEinfo.s);
	var res = parseFloat(HEinfo.res);
	var scents = 1200*log2(HEinfo.s+1);
	var padding = Math.round(100*scents);
	var min = parseFloat(HEinfo.mincents)-padding;
	var max = parseFloat(HEinfo.maxcents)+padding;
	var dist = HEinfo.dist;
	var series = HEinfo.series;
	var normalize = HEinfo.normalize;
	var locr = r;
	
	if(dist == 'linear') {
		min = 0-padding;
		max = 1200+padding;
	}

	a = (a == 1) ? a = 1.0000000001 : a = a;
	var rcount = 0;
	
	//build kernel
	var k = new Array();
	var ak = new Array();
	
	for(var i=(max-min)/res;i>=0;i--) {
		k[i] = 0;
		ak[i] = 0;
	}
	
	for(var i=locr.length-1; i>=0; i--) {
		var rcent, rcompl;
		
		if(dist == "log")
			rcent = 1200*log2(locr[i][0] / locr[i][1]);
		else if (dist == "linear")
			rcent = 1200*locr[i][0] / locr[i][1];
		
		if(series == "tenney")
			rcompl = Math.sqrt(locr[i][0] * locr[i][1]);
		else if (series == "farey")
			rcompl = locr[i][1];

		//check for bounds to optimize
		if(rcent < min || rcent > max)
			continue;
		
		rcount++;
		
		//start building kernel, first check for rounded off case that doesn't need interpolation
		if(rcent == Math.round(rcent)) {
			var index = (rcent - min)/res;
			k[index] += 1/rcompl;
			ak[index] += 1/Math.pow(rcompl,a);
		}
		//or else we do need interpolation
		else{
			var clow = Math.ceil(rcent) - rcent;
			var chigh = rcent - Math.floor(rcent);
			
			var index = (Math.floor(rcent) - min)/res;
			k[index] += 1/rcompl * clow;
			k[index+1] += 1/rcompl * chigh;
			
			ak[index] += 1/Math.pow(rcompl,a) * clow;
			ak[index+1] += 1/Math.pow(rcompl,a) * chigh;
		}
	}		


	//do convolution
	//first pad to a power of two to make the convolution easier for numerous reasons
	var oldlen = k.length;
	var minlen = 1;
	while (minlen	< 2*k.length)	//2 is to avoid circular convolution distortion
		minlen *= 2;

	for(var i=k.length;i<minlen;i++) {
		k[i] = 0;
		ak[i] = 0;
	}

	//now work out gaussian
	var g = new Array(minlen);
	var ag = new Array(minlen);
	var g_sum = 0;
	var ag_sum = 0;
	for(var i=minlen-1;i>=0;i--) {
		var c = i*res+min;
		var gval = (1/(scents*2*Math.PI)) * Math.exp(-(Math.pow(c-min,2)/(2*scents*scents))) + (1/(scents*2*Math.PI)) * Math.exp(-(Math.pow(c-(minlen*res+min),2)/(2*scents*scents)));
		g[i] = gval;
		g_sum += gval;
	}
	//normalize gaussian
	for(var i=g.length-1; i>=0; i--) {
		g[i] /= g_sum;
		ag[i] = Math.pow(g[i],a);
		ag_sum += ag[i];
	}

	var ent = conv(ak, ag);
	var nrm = conv(k, g);
	
	var nrmfct;
	if(normalize) {
		nrmfct = Math.log(rcount);
	}
	else
		nrmfct = 1;

	//trim answer and out
	var out = new Array();
	for (var i=oldlen-padding/res-1; i>=padding/res; i--) {
		var outval = ((1/(1-a)) * Math.log(ent[i]/Math.pow(nrm[i], a))/nrmfct);
		out[i-padding/res] = [(i*res + min), outval];
	}
	gg = g;
	gag = g;
	gk = k;
	gak = ak;
	gent = ent;
	gnrm = nrm;
	gout = out;
	return out;
}

function plotHE(HEinfo) {
	$(function () {
		var HE = harmonicEntropy(HEinfo);
		var tf = function(axis) {
			var t = new Array();
			for(i = axis.min; i-axis.max<=(axis.max-axis.min)/16; i+= (axis.max-axis.min)/8)
				t.push(i);
				
			if(t.length == 10)
				t.pop();
			return t;
		}
	  plot = $.plot($("#HEPlotDiv"), [{color:"#00C000", data: HE}], {xaxis: {tickSize: (HEinfo.maxcents-HEinfo.mincents)/12}, yaxis:{ticks:tf, autoscaleMargin:0.03, labelWidth:30, reserveSpace:true}, grid: {hoverable: true, clickable: true}});
	
	  function showTooltip(id, cls, color, item, extrastr) {
	  	if(!extrastr)
	  		extrastr = "";
	  	var x = item.pageX;
	  	var y = item.pageY;
	  	var dyad = item.datapoint[0].toFixed(3);
	  	var ent = item.datapoint[1].toFixed(3);
	    var t = $('<div id="'+id+'" class="'+cls+'">Dyad: ' + dyad + ' cents\n<br>HE: ' + ent + ' nats '+extrastr+'</div>').css( {
	    	position: 'absolute',
	      display: 'none',
	      top: y + 5,
	      left: x + 5,
	      border: '1px solid #fdd',
	      font: '8pt sans-serif',
	      padding: '2px',
	      'background-color': color,
	      cursor: 'move',
	      opacity: 0.80
	    });
	    t.appendTo("body").show().draggable().dblclick(function() {t.remove();plot.unhighlight(item.series, item.datapoint)});
	  }
	
	  var previousPoint = null;
	  $("#HEPlotDiv").bind("plothover", function (event, pos, item) {
	    if (item) {
	  		if (previousPoint != item.dataIndex) {
			    previousPoint = item.dataIndex;
			    
			    $("#tooltip").remove();
			    showTooltip("tooltip", "tooltip", "#ffeeee", item);
	     	}
	   	}
			else {
	    	$("#tooltip").remove();
	      previousPoint = null;            
	  	}
		});
	
		$("#HEPlotDiv").bind("plotclick", function (event, pos, item) {
	    if (item) {
				var x = item.datapoint[0].toFixed(3);
				var y = item.datapoint[1].toFixed(3);

				var name = "permatipx"+x.toString().replace(/\./gi, "p") +"y"+y.toString().replace(/\./gi, "p");
				var prevPoints = $("#"+name);
				if(prevPoints.length == 0) {
					showTooltip(name, "permatip", "#ffbbbb", item, "<br><small><i>Double-click to remove</i></small>");
					plot.highlight(item.series, item.datapoint);
	      }
	    }
	  });
	});
}

function preCalcRatios(HEinfo) {
	if(HEinfo.N != oldN || HEinfo.series != oldseries) {
		$("#progress-dialog").dialog({
			height: 100,
			width: 400,
			modal: true
		});
		r.length=0;
		
		oldN = HEinfo.N;
		oldseries = HEinfo.series;
		
		var n = HEinfo.N;
	setTimeout(function(){
    var start = +new Date();

			if(HEinfo.series == "tenney") {
        do {
					for(var i=1; i<=Math.floor(Math.sqrt(n)); i++) {
						if(n/i == Math.round(n/i) && gcd(i,n/i) == 1) {
							r.push([i, n/i]);		//does numerator on left, denominator on right
							if(n/i != i)
								r.push([n/i, i]);
						}
					}
        } while (--n >= 0 && (+new Date() - start < 50));
      }
      else if (HEinfo.series == "farey") {
        do {
					for(var i=0; i<=n; i++) {
						if(gcd(i,n) == 1) {
							r.push([n, i]);
							if(n != i)
								r.push([i, n]);
						}
					}
        } while (--n >= 0 && (+new Date() - start < 50));
      }

    if (n > 0){
					$("#progressbar").progressbar("option","value",100*(HEinfo.N-n)/HEinfo.N)
        setTimeout(arguments.callee, 25);
    } else {
					$("#progressbar").progressbar("option","value",100)
					$("#progress-dialog").dialog("close");
        setTimeout(function(){plotHE(HEinfo)},1);
    }
	}, 25);
	}
	else {
		plotHE(HEinfo);
	}
}

// Frontend
$(function() {
    $( "#slidera" ).slider({
        value:1,
        min: 0,
        max: 7,
        step: 0.1,
        slide: function( event, ui ) {
            $("#aval").text("a: " + ui.value.toFixed(1))
        }
    });
    
    $( "#sliders" ).slider({
        value:1,
        min: 0.3,
        max: 2,
        step:0.05,
        slide: function( event, ui ) {
            $("#sval").text("s: " + ui.value.toFixed(2) + "%")
        }
    });
    
    $( "#sliderres" ).slider({
        value:0,
        min: -2,
        max: 0,
        slide: function( event, ui ) {
            $("#resval").text("resolution: " + Math.pow(10, ui.value).toFixed(2) + " cents")
        }
    });
    
    $("#seriesradio").buttonset();
    
    $("#distradio").buttonset();

    $("#go").button().click(function(){
        plotted = true;
        if(plot !== undefined)
            plot.unhighlight();
        $(".permatip").remove();
        var HEinfo = {
            N: parseInt($("#textN").val()),
            s: parseFloat($("#sliders").slider("option", "value")/100),
            a: parseFloat($("#slidera").slider("option", "value")),
            series: $("input[name=seriesradiobuts]:checked").attr("id"),
            dist: $("input[name=distradiobuts]:checked").attr("id"),
            mincents: parseFloat($("#textmin").val()),
            maxcents: parseFloat($("#textmax").val()),
            res: Math.pow(10,parseFloat($("#sliderres").slider("option", "value"))),
            normalize: ($("#normcheck:checked").val() !== undefined)
        }
        preCalcRatios(HEinfo);
    });

    $("#csvbut").button().click(function(){
        if(!plotted)
            return;
        var str = "";
        var pd = plot.getData()[0].data;

        for (var i=0;i<pd.length;i++) {
            str += pd[i][0] + "," + pd[i][1] + "\n";
        }
            
        $("#csv").text(str);
        $("#dialog-modal").dialog({
            height: 400,
            width: 800,
            modal: true
        });
    });
    
    $("#progressbar").progressbar({
        value: 0
    });
});
