import os

class JSDiv:
    def __init__(self):
        pass
    def div(self):
        return "<div>override JSDiv.div()</div>"
    def page(self):
        return """<html><head><title>%(title)s</title>
<!--
<link href="/qap/js/jquery-ui/css/ui-lightness/jquery-ui-1.8.20.custom.css" rel=
"stylesheet">
<script type="text/javascript" src="/qap/js/jquery-ui/js/jquery.js"></script>
<script type="text/javascript" src="/qap/js/jquery-ui/js/jquery-ui.js"></script>
-->
</head>
<body>
%(divtent)s
</body></html>
""" % {"title":type(self), "divtent":self.div()}

class JSAce(JSDiv):
    code = None
    lnum = 0
    local_client = False
    filename = None
    basename = None
    dirname = None
    def __init__(self, code = None, 
                        lnum = 0, 
                        local_client = False, 
                        filename=None, 
                        target=None):
        JSDiv.__init__(self)
        self.code = code
        self.lnum = lnum
        self.local_client = local_client
        self.filename = filename
        self.target = target
        self.basename = os.path.basename(filename)
        self.dirname  = os.path.dirname(filename)
    def div(self,code = None):
        code = self.code
        if self.local_client:
            local_client_frag = '<input style="float:right" type="submit" value="SAVE"/>'
        else:
            local_client_frag = "<span style='color:dk_blue'>To be able to save you must be on localhost.<br/></span>"
        if code == None:
            code = """class Foo:
    prop = None
    def __init__(self):
    prop = "name"    
        """

        return """

        <div id="ace_div">
<style type="text/css" media="screen">
    #editor { 
        position: absolute;
        top: 30;
        right: 0;
        bottom: 0;
        left: 0;
    }
</style>
<div style="font-size:12px;position:absolute; top:0; bottom:30; left:0;right:0">
%(localClientFrag)s
Target: %(target)s from file: %(basename)s <br/>
in %(dirname)s<br/>

</div>
<div id="editor" style="width:600px">%(code)s</div>
    
<script src="http://d1n0x3qji82z53.cloudfront.net/src-min-noconflict/ace.js" type="text/javascript" charset="utf-8"></script>
<script>
        var editor = ace.edit("editor");
        editor.setTheme("ace/theme/monokai");
        editor.getSession().setMode("ace/mode/python");
        window.resizeTo(600,500);        
        window.setTimeout(function ()
            {
               editor.gotoLine(%(lnum)d,0,true);
               editor.scrollToLine(%(lnum)d-1);
               
            }, 100);
        
</script></div>
        """  % {"code":code, 
                "lnum":self.lnum, 
                "localClientFrag": local_client_frag,              
                "localClient":self.local_client,
                "filename":self.filename,
                "dirname":self.dirname,
                "basename":self.basename,
                "target":self.target
                }
class JSAccord(JSDiv):
    def div(self):
        return """
            <div class="gem_accord" style="display:none">
            <script type="text/javascript">
            function dict2accordian(accord, appTo, startOpen )
            {
                if (typeof(startOpen) == "undefined")
                {
                    var block = "block";
                    var subblock = "none"; 
                    var substartopen = false;
                }
                else
                {
                    if (startOpen)
                    {
                        var block = "block";
                        var subblock = block; 
                        var substartopen = startOpen
                    }                        
                    else
                    {
                        var block        = "none";
                        var subblock     = block; 
                        var substartopen = startOpen;

                    }
                }
                if (typeof(appTo) == "undefined")
                    appTo = null;
                // turns nested dict into collapsable hierarchical set of divs
                // distinguishes between typeof == object for sub-directories
                // and typeof == string which is an element
                            /*    
                            { category: name
                              subcats: [ {} {} {} ]
                              members: [ {} {} {} ]
                            }
                            */  
                var ndiv_html = '<div style="'
                                    + 'display:'+block+';'
                                    + 'border:solid black 1px; padding:2px;margin:1px;margin-left:1em"' 
                                    + 'class="recipies_'
                                    + accord.cat_name
                                    + '">'
                                    + '</div>'
                var title_ndiv_html = '<div style="margin-left:1em"><a href="javascript:void(0)" onclick="'
                                + "$('.recipies_"
                                    +accord.cat_name
                                    +"').toggle()"
                                + '">'
                                + "Recipe Folder: "
                                + accord.cat_name
                                + '</a></div>'
                
                var title = $(title_ndiv_html);
                if (appTo)
                {
                    appTo.append(title);
                }                  
                //console.log("js25:"+ndiv_html);
                var ndiv = $(ndiv_html);
                //for (key in accord)
                //{
                //    console.log("js27:"+ key);
                //}
                if ("members" in accord)
                {
                    accord.members.sort( function (a,b) {
                        var aname  = a.name.toLowerCase();
                        var bname  = b.name.toLowerCase();
                        if (aname > bname) { return 1;};
                        if (aname < bname) { return -1;};                        
                        return 0;                        
                        });
                    for (var i = 0; i < accord.members.length; i++)
                    {
                        var mem = accord.members[i];
                        // console.log("js26:"+JSON.stringify(mem));
                        var mem_span ='<span title="'
                                            + mem.path
                                            + '"'
                                            + ' class="recipe_'
                                            + mem.name
                                            + '">'
                                            + ' (<a style = "font-size:75%" href="javascript: void(0)" onclick="editFile('
                                            + "'"
                                            + mem.path
                                            + "')"+ '">'
                                            + "visit</a>) "
                                            + mem.name 
                                            
                                            + '</span></br>'; 
                        //console.log("js69:"+mem_span);
                        ndiv.append($(mem_span));
                    }
                }
                if ("subcats" in accord)
                {
                    for (var i = 0; i < accord.subcats.length; i++)
                    {
                        var subcat = accord.subcats[i];
                        var nsdiv = dict2accordian(subcat, ndiv, substartopen);
                        ndiv.append(nsdiv);
                    }
                }
                //console.log("js28:"+ndiv.html());
                return ndiv;
            }
            </script>
            </div>
            """

class JSTypes(JSDiv):
    def div(self):
        return """<div class="lp_types" style="width:20%;float:left">
                   <script type="text/javascript">
                   function populateTypesDiv(typemap)
                   {
                        // console.log("ju10:"+JSON.stringify(typemap))
                        var typl = typemap["types"];
                        buff = "<b>Types within <tt>"
                             + typemap["package"]
                             + "</tt> Package</b><br/>";
                        for (var i =0 ; i < typl.length; i++)
                        {
                            buff += '(<a style ="font-size:75%"href="javascript: void(0)" onclick="editFile(';
                            buff += "'"
                            buff += typemap["type_meta"][typl[i]].fullpath;
                            buff += "')"+ '">';
                            buff += "visit"
                            buff += "</a>) ";
                            buff += typl[i]; 
                            /* server editor
                            buff += ' (<a href="javascript: void(0)" onclick="ajaxLink(';
                            buff += "'"
                            buff += typemap["type_meta"][typl[i]].edit_link;
                            buff += "')"+ '">';
                            buff += "visit"
                            buff += "</a>)";
                            */
                            buff += "<br/>";
                        }
                        //console.log("16:"+buff);
                        $($(".lp_types_content")[0]).html(buff);
                   }
                   </script>
                   <div class = "lp_types_content" style="border:solid black 1px;padding:2px;margin:3px">
                   </div>                               
                </div> 
                """
class BROKENJSDescriptors(JSDiv):
    def div(self):
         return """<div class="lp_descriptors" style="float:left;width:38%">
                   <script type="text/javascript">
                   function populateDescriptorsDiv(desmap)
                   {
                        console.log("js256:",desmap);
                        var descmap = desmap["descriptors"]
                        buff = "<b>Descriptor Calculators within <tt>"
                                +desmap["package"]
                                +"</tt> Package</b><br/>";
                        buff += "<div  style='overflow:hidden;width:38%;border:solid black 1px; margin:2px;padding:1px'>";
                        buff += "<div style='width:25%'>Type</div>"
                        buff += "<div style='float:left;'>Calculator Class</div>"
                        for (var i = 0; i < descmap.length; i++)
                        {
                            var desc = descmap[i];
                            
                            line =      "<div style='width:25%'>" + desc.type + "</div>"
                                        +"<div style='float:left;'><span title='"
                                        + desc.path
                                        + "'>" 
                                        + desc.descriptor_class
                                        + "</span></div>";
                            //console.log("47:"+line);
                            buff+=line;
                        }
                        buff += "</div>";
                        $($(".lp_descriptors_content")[0]).html(buff);
                        
                   }
                   </script>
                   <div class = "lp_descriptors_content" style="border:solid black 1px;padding:2px;margin:3px">
                   </div>                              
                </div> 
                """     
   
class JSDescriptors(JSDiv):
    def div(self):
         return """<div class="lp_descriptors" style="float:left;width:38%">
                   <script type="text/javascript">
                   function showCalculator(astrotype)
                   {
                        var loaded = false;
                        var focus = $($(".pdk_focus")[0]);
                        focus.slideUp(500,function () {
                            if (!loaded){
                                focus.empty();
                                focus.html("Loading Calculator Information: "  + astrotype);                            
                                function pulse ()
                                {
                                    if (!loaded)
                                    {
                                        console.log("pulse");
                                        focus.fadeOut(500).fadeIn(500);
                                        setTimeout(pulse,1000);
                                    }                                
                                }
                                focus.fadeIn(1000)
                                setTimeout( pulse, 1000);
                                    
                                }                            
                            });
                        $.ajax({url: "/calculator_by_type/"+astrotype,
                                data:{  "astrotype":astrotype},
                                type:"get",
                                dataType:"json",
                                success: function(data) {
                                    loaded = true;
                                    focus.empty();
                                    focus.slideUp(500, function () {
                                        
                                        
                                        var jqf = 'empty_pdk_focus()';
                                        var html =  "<span>Calculator For Type: <b><i>"
                                                    + data.astrotype 
                                                    + "</i></b> "
                                                    + "<a style = 'font-size:70%%' href='javascript:void(0)' onclick='" 
                                                    + jqf
                                                    + "'>(clear)</a>"
                                                    
                                                    + "</span>";
                                        var descriptors = data.descriptorDict;
                                        
                                        focus.append($(html));
                                        if (descriptors.length == 0)
                                        {
                                        focus.append("NO Descriptors DEFINED");
                                        focus.slideDown();
                                        }
                                        else
                                        {
                                            //console.log("215:"+JSON.stringify(data.prim_info));
                                            for ( key in descriptors)
                                            {
                                                
                                                var lnum  = descriptors[key]["lnum"];
                                                var tpath = descriptors[key]["path"];
                                                var name = key
                                                         
                                                primml = "<div style='margin-left:5em'>"
                                                         + " (<a style='font-size:75%' href='javascript:void(0)' "
                                                         + " onclick='aceditFile("
                                                         + '"' + tpath + '", '
                                                         + '"' + lnum + '")'
                                                         + "'>"
                                                         + "visit"
                                                         + "</a>) "
                                                         + "<span title='"+tpath+"'>"
                                                         + key
                                                         + "</span>"
                                                         + "</div>";
                                                focus.append($(primml));
                                            } 
                                            focus.slideDown(); 
                                            $("html, body").animate({ scrollTop: 0 }, 1000);               
                                        }                                    
                                    }); // slideup func
                                 }  // ajax dict                             
                                }); //top slideup                    
                   }
                   
                   function populateDescriptorsDiv(desmap)
                   {
                        //alert(JSON.stringify(desmap));
                        var descmap = desmap["descriptors"]
                        buff = "<b>Descriptor Calculators within <tt>"
                                +desmap["package"]
                                +"</tt> Package</b><br/>";
                        buff += "<table  width = '38%' border='1px' cellspacing='2px' cellpadding='1px'>";
                        buff += "<tr><th>Type</th><th>Calculator Class</th></tr>"
                        for (var i = 0; i < descmap.length; i++)
                        {
                            var desc = descmap[i];
                            
                            line =      "<tr>"
                                        +"<td>" + desc.type + "</td>"
                                        +"<td><span title='"
                                        + desc.path
                                        + "'>" 
                                        + "("
                                        + '<a style="font-size:75%" href="javascript:void(0)" onclick="showCalculator('
                                        + "'" + desc.type + "')" + '">'
                                        + "detail"
                                        + "</a>"
                                        + ") "
                                        + desc.descriptor_class
                                        + "</span></td></tr>";
                            //console.log("47:"+line);
                            buff+=line;
                        }
                        buff += "</table>";
                        var targel = $($(".lp_descriptors_content")[0])
                        targel.css("overflow", "hidden")                        
                        targel.html(buff);
                        
                   }
                   </script>
                   <div class = "lp_descriptors_content" style="border:solid black 1px;padding:2px;margin:3px">
                   </div>                              
                </div> 
                """     
   
class JSRecipeSystem(JSDiv):
    def div(self):
         return """<div class="lp_recipesystem" style="float:left;width:38%">
                   <script type="text/javascript">
                   function showPrimset(primclass, astrotype)
                   {
                        var loaded = false;
                        focus = $($(".pdk_focus")[0])
                        focus.slideUp(500,function () {
                            if (!loaded){
                                focus.empty();
                                focus.html("Loading Primset Information: "  + primclass);                            
                                function pulse ()
                                {
                                    if (!loaded)
                                    {
                                        console.log("pulse");
                                        focus.fadeOut(500).fadeIn(500);
                                        setTimeout(pulse,1000);
                                    }                                
                                }
                                focus.fadeIn(1000)
                                setTimeout( pulse, 1000);
                                    
                                }                            
                            });
                        $.ajax({url: "/primset_by_type/"+astrotype,
                                data:{ primclass: primclass,
                                        "astrotype":astrotype},
                                type:"get",
                                dataType:"json",
                                success: function(data) {
                                    loaded = true;
                                    focus.empty();
                                    focus.slideUp(500, function () {
                                        
                                        var otherfrag = "";
                                        if (data.other_primsets.length)
                                        {
                                            for (var i = 0; i<data.other_primsets.length; i++)
                                            {
                                                otherfrag += data.other_primsets[i].link_frag;
                                                otherfrag += " ";
                                            }
                                            otherfrag = "("+otherfrag+")";
                                        }
                                        
                                        var jqf = 'empty_pdk_focus()';
                                        var html =  "<span>Primitive Set Class: <b><i>"
                                                    + data.class 
                                                    + "</i></b> "
                                                    + "<a style = 'font-size:70%%' href='javascript:void(0)' onclick='" 
                                                    + jqf
                                                    + "'>(clear)</a>"
                                                    
                                                    + "<span style='font-size:75%'> "
                                                    + otherfrag
                                                    + "</span>"
    
                                                    + "<br/><u>from <tt>"
                                                    + data.path
                                                    + "</tt><br/>"
                                                    + "</u></span>";
                                        var prims = data.prims;
                                        var pinfo = data.prim_info;
                                        
                                        
                                        focus.append($(html));
                                        if (prims.length == 0)
                                        {
                                        focus.append("NO PRIMITIVES DEFINED");
                                        focus.slideDown();
                                        }
                                        else
                                        {
                                            //console.log("215:"+JSON.stringify(data.prim_info));
                                            var basepath = data["path"];
                                            for (var i =0; i < prims.length; i++)
                                            {
                                                var lnum  = data.prim_info[prims[i]]["lnum"];
                                                var tpath = data.prim_info[prims[i]]["path"];
                                                // console.log("382:"+tpath+"<->"+basepath);
                                                if (basepath == tpath)
                                                {
                                                    fromstr = "";
                                                }
                                                else
                                                {
                                                    fromstr = "<tt style='font-size:70%;color:gray' "
                                                                + "title='"
                                                                + data.prim_info[prims[i]]["basename"]
                                                                + "'> inherited "
                                                                +"</tt>";
                                                }                      
                                                primml = "<div style='margin-left:5em'>"
                                                         + " (<a style='font-size:75%' href='javascript:void(0)' "
                                                         + " onclick='aceditFile("
                                                         + '"' + tpath + '", '
                                                         + '"' + lnum + '")'
                                                         + "'>"
                                                         + "visit"
                                                         + "</a>) "
                                                         + "<span title='"+tpath+"'>"
                                                         + prims[i]
                                                         + "</span>"
                                                         + fromstr
                                                         + "</div>";
                                                focus.append($(primml));
                                            } 
                                            focus.slideDown(); 
                                            $("html, body").animate({ scrollTop: 0 }, 1000);               
                                        }                                    
                                    }); // slideup func
                                 }  // ajax dict                             
                                }); //top slideup                    
                   }
                   function populateRecipeSystemDiv(fullmap)
                   {
                        //alert(JSON.stringify(recmap));
                        if (false)
                        {
                            var     recmap = fullmap["recipes"]
                            buff = "<b>Recipies within <tt>"
                                    +fullmap["package"]
                                    +"</tt> Package</b><br/>";
                            for (var recname in recmap)
                            {
                                var rec = recmap[recname];
                                //console.log("js86:"+rec)
                                line =      "<span title='"
                                            + rec.path
                                            + "'>" 
                                            + rec.name 
                                            + "</span><br/>";
                                // console.log("47:"+line);
                                buff+=line;
                            }
                            $($(".lp_recipes_content")[0]).html(buff);
                        }
                        var el = $($(".lp_recipes_content")[0]);
                        el.empty();
                        el.append("<b>Recipes Declared within <tt>"
                                    +fullmap["package"]
                                    +"</tt> Package</b><br/>");
                        el.append(
                                dict2accordian(fullmap["recipe_cat"], el)
                                );
                        // console.log("178:"+JSON.stringify(fullmap["recipe_cat"], null, 3));
                        
                        
                        // PRIMITIVE_SET
                        var prim_cat = fullmap["primitives_cat"];
                        el = $($(".lp_primitives_content")[0]);
                        el.empty();                        
                        var keys = getKeys(prim_cat);
                        
                        buff = "<table border='1px' cellspacing='2px' cellpadding='1px' >";
                        buff += "<tr><th>Type</th><th>Primitive Set Class</th></tr>"
                        for (var i = 0; i < keys.length; i++)
                        {
                            var key = keys[i];
                            var ps = prim_cat[key];
                            line =      "<tr>"
                                        +"<td>" + key + "</td>"
                                        +"<td><span title='"
                                        + "from index: "  + ps.index_path
                                        + "'>" 
                                        + "(<a style='font-size:75%' href='javascript:void(0)'"
                        	            + " onclick='showPrimset("
                        	            + '"' + ps.class + '","'
                        	            + ps.astrotype + '")' + "'>"
                        	            + "detail</a>) "
                        	            
                        	            + ps.class
                        	                                                  
                                        + "</span></td></tr>";
                            //console.log("47:"+line);
                            buff+=line;
                        }
                        buff += "</table>";
                        el.append($(buff));
                   }
                   </script>
                   <div class = "lp_recipes_content" style="border:solid black 1px;padding:2px;margin:3px">
                   </div>
                   <div class = "lp_primitives_content" style="border:solid black 1px;padding:2px;margin:3px">
                   </div>                               
                </div> 
                """      
class JSPackageSelection(JSDiv):
    options = None
        
    def __init__(self, opts):
        self.options = opts

    def div(self):
        ret = """
            <script type="text/javascript">
            function loadPkgContent(url, name)
            {
                console.log("js613: here", url);
                $(".pdk_focus").slideUp();
                $.ajax({ url: url,
                        dataType: "json",
                        success: function(data, ts, jq)
                                 {
                                     console.log(JSON.stringify(data));
                                     populateTypesDiv(data);
                                     populateDescriptorsDiv(data);
                                     populateRecipeSystemDiv(data);
                                 }
                        }
                    );
            }            
            </script>
            """ 
        # ret += "<ul>"
        
        for el in self.options:
            ret += """<input type="submit" value="%(name)s" onclick="loadPkgContent('%(url)s', '%(name)s')"></input>
                      <br/>
                   """ % {
                          "name" : el["name"], 
                          "url"  : el['url']
                         }
            
        tret = """  <div class="adk_packages" style="padding:5px;border:solid black 1px;float:right">
                    <u>packages</u><br/>
                    %s                    
                    </div>              
                """ % ret
        return tret
