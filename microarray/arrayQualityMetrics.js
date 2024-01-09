// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, true, false, false, false, false, false, false, true, false, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, true, false, true, false, false, false, false ];
var arrayMetadata    = [ [ "1", "164_I_.CEL", "164_I", "164", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "2", "164_II.CEL", "164_II", "164", "Crohn's disease", "inflamed colonic mucosa" ], [ "3", "183_I.CEL", "183_I", "183", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "4", "183_II.CEL", "183_II", "183", "Crohn's disease", "inflamed colonic mucosa" ], [ "5", "2114_I.CEL", "2114_I", "2114", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "6", "2114_II.CEL", "2114_II", "2114", "Crohn's disease", "inflamed colonic mucosa" ], [ "7", "2209_A.CEL", "2209_A", "2209", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "8", "2209_B.CEL", "2209_B", "2209", "Crohn's disease", "inflamed colonic mucosa" ], [ "9", "2255_I.CEL", "2255_I", "2255", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "10", "2255_II.CEL", "2255_II", "2255", "Crohn's disease", "inflamed colonic mucosa" ], [ "11", "2400_I.CEL", "2400_I", "2400", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "12", "2400_II.CEL", "2400_II", "2400", "ulcerative colitis", "inflamed colonic mucosa" ], [ "13", "2424_A.CEL", "2424_A", "2424", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "14", "2424_B.CEL", "2424_B", "2424", "ulcerative colitis", "inflamed colonic mucosa" ], [ "15", "255_I.CEL", "255_I", "255", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "16", "255_II.CEL", "255_II", "255", "Crohn's disease", "inflamed colonic mucosa" ], [ "17", "2826_I.CEL", "2826_I", "2826", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "18", "2826_II.CEL", "2826_II", "2826", "Crohn's disease", "inflamed colonic mucosa" ], [ "19", "2853_I.CEL", "2853_I", "2853", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "20", "2853_II.CEL", "2853_II", "2853", "Crohn's disease", "inflamed colonic mucosa" ], [ "21", "2978_I.CEL", "2978_I", "2978", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "22", "2978_II.CEL", "2978_II", "2978", "Crohn's disease", "inflamed colonic mucosa" ], [ "23", "2987_I.CEL", "2987_I", "2987", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "24", "2987_II.CEL", "2987_II", "2987", "ulcerative colitis", "inflamed colonic mucosa" ], [ "25", "2992_I.CEL", "2992_I", "2992", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "26", "2992_II.CEL", "2992_II", "2992", "ulcerative colitis", "inflamed colonic mucosa" ], [ "27", "2995_I.CEL", "2995_I", "2995", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "28", "2995_II.CEL", "2995_II", "2995", "ulcerative colitis", "inflamed colonic mucosa" ], [ "29", "321_I.CEL", "321_I", "321", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "30", "321_II.CEL", "321_II", "321", "Crohn's disease", "inflamed colonic mucosa" ], [ "31", "3222_I.CEL", "3222_I", "3222", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "32", "3222_II.CEL", "3222_II", "3222", "ulcerative colitis", "inflamed colonic mucosa" ], [ "33", "3223_I.CEL", "3223_I", "3223", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "34", "3223_II.CEL", "3223_II", "3223", "ulcerative colitis", "inflamed colonic mucosa" ], [ "35", "3226_I.CEL", "3226_I", "3226", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "36", "3226_II.CEL", "3226_II", "3226", "ulcerative colitis", "inflamed colonic mucosa" ], [ "37", "3233_I.CEL", "3233_I", "3233", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "38", "3233_II.CEL", "3233_II", "3233", "ulcerative colitis", "inflamed colonic mucosa" ], [ "39", "3258_I.CEL", "3258_I", "3258", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "40", "3258_II.CEL", "3258_II", "3258", "ulcerative colitis", "inflamed colonic mucosa" ], [ "41", "3259_I.CEL", "3259_I", "3259", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "42", "3259_II.CEL", "3259_II", "3259", "ulcerative colitis", "inflamed colonic mucosa" ], [ "43", "3262_I.CEL", "3262_I", "3262", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "44", "3262_II.CEL", "3262_II", "3262", "Crohn's disease", "inflamed colonic mucosa" ], [ "45", "3266_I.CEL", "3266_I", "3266", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "46", "3266_II.CEL", "3266_II", "3266", "Crohn's disease", "inflamed colonic mucosa" ], [ "47", "3269_I.CEL", "3269_I", "3269", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "48", "3269_II.CEL", "3269_II", "3269", "ulcerative colitis", "inflamed colonic mucosa" ], [ "49", "3271_I.CEL", "3271_I", "3271", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "50", "3271_II.CEL", "3271_II", "3271", "Crohn's disease", "inflamed colonic mucosa" ], [ "51", "3302_I.CEL", "3302_I", "3302", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "52", "3302_II.CEL", "3302_II", "3302", "Crohn's disease", "inflamed colonic mucosa" ], [ "53", "3332_I.CEL", "3332_I", "3332", "Crohn's disease", "non-inflamed colonic mucosa" ], [ "54", "3332_II.CEL", "3332_II", "3332", "Crohn's disease", "inflamed colonic mucosa" ], [ "55", "848_A.CEL", "848_A", "848", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "56", "848_B.CEL", "848_B", "848", "ulcerative colitis", "inflamed colonic mucosa" ], [ "57", "888_I.CEL", "888_I", "888", "ulcerative colitis", "non-inflamed colonic mucosa" ], [ "58", "888_II.CEL", "888_II", "888", "ulcerative colitis", "inflamed colonic mucosa" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
