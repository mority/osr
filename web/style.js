// Map style: the MOTIS UI style (Shortbread vector tiles), see motis-style.js.
// Tiles, glyphs and sprites all come from europe.motis-project.de (CORS
// enabled for all of them; Transitous lacks CORS on glyphs and sprites).
import {getStyle} from "./motis-style.js";

const kBaseUrl = 'https://europe.motis-project.de/';

const darkQuery = window.matchMedia('(prefers-color-scheme: dark)');

function isDark() {
    return darkQuery.matches;
}

function onThemeChange(fn) {
    darkQuery.addEventListener('change', fn);
}

function style(map, level) {
    map.setStyle(getStyle(isDark() ? 'dark' : 'light', level, kBaseUrl, kBaseUrl,
        false));
}

export {style, isDark, onThemeChange}
