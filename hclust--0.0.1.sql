CREATE OR REPLACE FUNCTION
hclust(text,integer, text, integer, integer, text) RETURNS integer AS 'MODULE_PATHNAME','hclust'
LANGUAGE C STRICT;


