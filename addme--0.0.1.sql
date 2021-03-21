CREATE OR REPLACE FUNCTION
addme(text,integer, text, integer, integer, text) RETURNS integer AS 'MODULE_PATHNAME','addme'
LANGUAGE C STRICT;


