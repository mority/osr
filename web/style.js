import {getStyle} from "./motis-style.js";

// Tiles come from a MOTIS instance. The style is MOTIS's own
// (motis-style.js), because it has to match the tile schema that instance
// serves.
const kMotisBaseUrl = "https://api.transitous.org/";

// Glyphs and sprites are served from here: MOTIS sends CORS headers on its
// API only, so a browser refuses its static files on any other origin. They
// are copies of motis ui/static - glyphs/ is too big to commit and ignored.
const kStaticBaseUrl = `${window.location.origin}/`;

function createShield(opt) {
    const d = 32

    const cv = document.createElement('canvas');
    cv.width = d;
    cv.height = d;
    const ctx = cv.getContext('2d');

    // coord of the line (front = near zero, back = opposite)
    const l_front = 1;
    const l_back = d - 1;

    // coord start of the arc
    const lr_front = l_front + 2;
    const lr_back = l_back - 2;

    // control point of the arc
    const lp_front = l_front + 1;
    const lp_back = l_back - 1;

    let p = new Path2D();
    p.moveTo(lr_front, l_front);

    // top line
    p.lineTo(lr_back, l_front);
    // top right corner
    p.bezierCurveTo(lp_back, lp_front, lp_back, lp_front, l_back, lr_front);
    // right line
    p.lineTo(l_back, lr_back);
    // bottom right corner
    p.bezierCurveTo(lp_back, lp_back, lp_back, lp_back, lr_back, l_back);
    // bottom line
    p.lineTo(lr_front, l_back);
    // bottom left corner
    p.bezierCurveTo(lp_front, lp_back, lp_front, lp_back, l_front, lr_back);
    // left line
    p.lineTo(l_front, lr_front);
    // top left corner
    p.bezierCurveTo(lp_front, lp_front, lp_front, lp_front, lr_front, l_front);

    p.closePath();

    ctx.fillStyle = opt.fill;
    ctx.fill(p);
    ctx.strokeStyle = opt.stroke;
    ctx.stroke(p);

    return [
        ctx.getImageData(0, 0, d, d),
        {
            content: [lr_front, lr_front, lr_back, lr_back],
            stretchX: [[lr_front, lr_back]],
            stretchY: [[lr_front, lr_back]]
        }
    ];
}

// The level is baked into the style's filters, so a level change sets the
// style again - which drops every source and layer added on top of it.
const style = (map, level) => {
    map.setStyle(getStyle("light", level, kStaticBaseUrl, kMotisBaseUrl, false));

    // Road shields are drawn at runtime, as the MOTIS UI does.
    const shields = [
        ["shield", "hsl(0, 0%, 98%)", "hsl(0, 0%, 75%)"],
        ["shield-dark", "hsl(0, 0%, 16%)", "hsl(0, 0%, 30%)"]
    ];
    for (const [id, fill, stroke] of shields) {
        if (!map.hasImage(id)) {
            map.addImage(id, ...createShield({fill, stroke}));
        }
    }
};

export {style}
