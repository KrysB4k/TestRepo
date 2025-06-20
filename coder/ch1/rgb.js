"use strict";
let canvas = document.getElementById("rgb");
canvas.width = 256;
canvas.height = 256;

let ctx = canvas.getContext("2d");
ctx.fillStyle = "black";
ctx.fillRect(0, 0, canvas.width, canvas.height);

function setColor(red, green, blue)
{
  const color = `rgb(${red}, ${green}, ${blue})`;
  ctx.fillStyle = color;
  ctx.fillRect(0, 0, canvas.width, canvas.height);
}

function getRed()
{
  const slide = document.getElementById("slider-r");
  let r = slide.value;

  const txt = document.getElementById("val-r");
  txt.innerHTML = r.toString();

  return r;
}

function getGreen()
{
  const slide = document.getElementById("slider-g");
  let g = slide.value;

  const txt = document.getElementById("val-g");
  txt.innerHTML = g.toString();

  return g;
}

function getBlue()
{
  const slide = document.getElementById("slider-b");
  let b = slide.value;

  const txt = document.getElementById("val-b");
  txt.innerHTML = b.toString();

  return b;
}

function manageSliders()
{
  setColor(getRed(), getGreen(), getBlue());
}

manageSliders();

let reds = document.getElementById("slider-r");
reds.addEventListener("change", manageSliders);
reds.addEventListener("input", manageSliders);

let greens = document.getElementById("slider-g");
greens.addEventListener("change", manageSliders);
greens.addEventListener("input", manageSliders);

let blues = document.getElementById("slider-b");
blues.addEventListener("change", manageSliders);
blues.addEventListener("input", manageSliders);
