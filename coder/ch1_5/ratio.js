"use strict";
let canv = document.getElementById("ratio");
canv.width = 320;
canv.height = 320;

let ctx2 = canv.getContext("2d");
ctx2.fillStyle = "white";
ctx2.fillRect(0, 0, canv.width, canv.height);

function drawRect(xratio, yratio)
{
  ctx2.fillStyle = "white";
  ctx2.fillRect(0, 0, canv.width, canv.height);

  const xpoint = canv.width / 2;
  const ypoint = canv.height / 2;
  const width = xpoint * xratio * 2;
  const height = xpoint * yratio * 2;

  ctx2.fillStyle = "black";
  ctx2.fillRect(xpoint - width * 0.5, ypoint - height * 0.5, width, height);
}

function getRatio()
{
  const slide = document.getElementById("rat");
  let r = slide.value / 100;

  const txt = document.getElementById("val-rat");
  txt.innerHTML = r.toString();

  return r;
}

function updateSizeXY()
{
  const r = getRatio();
  const xratio = (r + 1) / 2;
  const yratio = 1 - xratio;

  const txt_x = document.getElementById("val-x");
  txt_x.innerHTML = (Math.round(xratio * 100) / 100).toString();

  const txt_y = document.getElementById("val-y");
  txt_y.innerHTML = (Math.round(yratio * 100) / 100).toString();

  drawRect(xratio, yratio);
}

updateSizeXY();

let rat = document.getElementById("rat");
rat.addEventListener("change", updateSizeXY);
rat.addEventListener("input", updateSizeXY);
