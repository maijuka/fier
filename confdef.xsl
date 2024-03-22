<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:template match="root">
<html> 
<head>
<title>SIREN &#8212; Interactive and visual geospatial redescription mining</title>
<style type="text/css">
body
{
margin:10px;
font-family:verdana,helvetica,sans-serif;
}

h3
{
display:inline;
}

div.page-header
{
padding:10px;
border:1px dotted #cccccc;
background-color:#efefef;
font-style:italic;
}

.parameter-name
{
margin-left:5px;
}

.parameter-legend
{
display:block;
color:#555555;
font-style:italic;
width:75%;
}
</style>
</head>
<body>
<h2>Preferences parameters</h2>

<div class="page-header">
<p><xsl:value-of select="info"/></p>

<p>
This is a list of available parameters listed by section in the form: 
<span class="parameter-name">[parameter name = default value]</span>
</p>
<p>
Parameters can be set via the interface menu Edit &#8594; Preferences.
</p>
<p>
An XML preferences file contains all non-default parameters. The ouline is as follows:
<pre class="xml-example">
&#60;root&#62;
	&#60;parameter&#62;
		&#60;name&#62;parameter name&#60;/name&#62;
		&#60;value&#62;parameter value&#60;/value&#62;
	&#60;/parameter&#62;
	&#60;parameter&#62;
		&#60;name&#62;parameter name&#60;/name&#62; 
		&#60;value&#62;parameter value&#60;/value&#62;
		&#60;value&#62;parameter value&#60;/value&#62; --- multiple choices parameters might have several values
	&#60;/parameter&#62;
	...
&#60;/root&#62;
</pre>
</p>

</div>

<xsl:for-each select="section">
	<h2 class="section"><xsl:value-of select="name"/></h2>
	<ul>
	<xsl:for-each select="section/parameter[parameter_type='open']">
		<li>  		
		<h3 class="parameter-label"><xsl:value-of select="label"/></h3>
		<span class="parameter-name">[<xsl:value-of select="name"/> =
		<xsl:value-of select="default/value"/>]
		</span>
		<p class="parameter-legend"><xsl:value-of select="legend"/>
		(<span class="parameter-details"><xsl:value-of select="value_type"/></span>)
		</p>
		</li>
	</xsl:for-each>
	<xsl:for-each select="section/parameter[parameter_type='range']">
		<li>  		
		<h3 class="parameter-label"><xsl:value-of select="label"/></h3>
		<span class="parameter-name">[<xsl:value-of select="name"/> =
		<xsl:value-of select="default/value"/>]
		</span>
		<p class="parameter-legend"><xsl:value-of select="legend"/>
		(<span class="parameter-details"><xsl:value-of select="value_type"/> in range [<xsl:value-of select="range_min"/> , <xsl:value-of select="range_max"/>]</span>)
		</p>		
		</li>
	</xsl:for-each>
        <xsl:for-each select="section/parameter[parameter_type='boolean']">
		<li>  		
		<h3 class="parameter-label"><xsl:value-of select="label"/></h3>
		<span class="parameter-name">[<xsl:value-of select="name"/> =
		<xsl:value-of select="default/value"/>]
		</span>
		<p class="parameter-legend"><xsl:value-of select="legend"/>
		(<span class="parameter-details">Boolean, yes/no</span>)
		</p>
		</li>
	</xsl:for-each>
	<xsl:for-each select="section/parameter[parameter_type='single_options']">
		<li>  		
		<h3 class="parameter-label"><xsl:value-of select="label"/></h3>
		<span class="parameter-name">[<xsl:value-of select="name"/> =
		<xsl:value-of select="default/value"/>]
		</span>
		<p class="parameter-legend"><xsl:value-of select="legend"/>
		(<span class="parameter-details">Single choice among: 
		<xsl:for-each select="options/value">
			<xsl:value-of select="."/>
			<xsl:if test="position() != last()">, </xsl:if>
		</xsl:for-each>
		</span>)
		</p>
		</li>
	</xsl:for-each>
	<xsl:for-each select="section/parameter[parameter_type='multiple_options']">
		<li>  		
		<h3 class="parameter-label"><xsl:value-of select="label"/></h3>
		<span class="parameter-name">[<xsl:value-of select="name"/> =
		<xsl:for-each select="default/value">
			<xsl:value-of select="."/>
			<xsl:if test="position() != last()">, </xsl:if>
		</xsl:for-each>]
		</span>
		<p class="parameter-legend"><xsl:value-of select="legend"/>
		(<span class="parameter-details">Multiple choices among: 
		<xsl:for-each select="options/value">
			<xsl:value-of select="."/>
			<xsl:if test="position() != last()">, </xsl:if>
		</xsl:for-each>
		</span>)
		</p>
		</li>
	</xsl:for-each>
	</ul>
</xsl:for-each>

</body>
</html>
</xsl:template>

</xsl:stylesheet>
