#!/bin/sh

headers="
lib/preamble.lua
lib/pixel.lua
lib/turtle.lua
templates/pixel_program.lua
templates/turtle_program.lua
examples/arbre.lua
examples/flake.lua
examples/sierpinski.lua
examples/sponge.lua
examples/creeper.lua
examples/plotter.lua
examples/maison.lua
games/hackman.lua
games/labyrinthe.lua
games/asteroids.lua
book/S01E01.lua
book/S01E02.lua
"

cat <<EOF
/*
 * This file was automatically generated, do not edit.
 */

#include <geogram/lua/lua_io.h>

void register_embedded_lua_files(void);
   
void register_embedded_lua_files() {
EOF

for header in $headers
do
    echo "     register_embedded_lua_file(\"$header\","
    cat $header | sed -e 's|\\n|\\\\n|' \
                      -e 's|"|\\\"|g' \
		      -e 's|^|        \"|' \
		      -e 's| *$| \\n\"|' 
    echo "     );"
    echo
done

echo "}"
