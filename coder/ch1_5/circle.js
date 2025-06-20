"use strict";
let canvas = document.getElementById("circle");
canvas.width = 320;
canvas.height = 320;

let ctx = canvas.getContext("2d");
ctx.fillStyle = "white";
ctx.fillRect(0, 0, canvas.width, canvas.height);

function drawArc(radius, angle)
{
  ctx.fillStyle = "white";
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  const a = angle * Math.PI / 180;
  ctx.beginPath();
  ctx.fillStyle = "blue";
  ctx.strokeStyle = "blue";
  ctx.lineWidth = 5;
  ctx.arc(160, 160, radius, 0, a, false);
  ctx.stroke();
  ctx.closePath();
  ctx.beginPath();
  ctx.fillStyle = "black";
  ctx.strokeStyle = "black";
  ctx.arc(160, 160, 6, 0, 2 * Math.PI);
  ctx.fill();
  ctx.closePath();
  ctx.beginPath();
  ctx.moveTo(160, 160);
  ctx.lineTo(160 + radius * Math.cos(a), 160 + radius * Math.sin(a));
  ctx.lineWidth = 1.5;
  ctx.stroke();
  ctx.closePath();
}

function getRadius()
{
  const slide = document.getElementById("radius");
  let r = slide.value;

  const txt = document.getElementById("val-r");
  txt.innerHTML = (r / 20).toString();

  return r;
}

function getAngle()
{
  const slide = document.getElementById("angle");
  let a = slide.value;

  const txt = document.getElementById("val-ang");
  txt.innerHTML = a.toString();

  return a;
}

function manageSliders()
{
  const radius = getRadius();
  const angle = getAngle();
  drawArc(radius, angle);
}

manageSliders();

let rad = document.getElementById("radius");
rad.addEventListener("change", manageSliders);
rad.addEventListener("input", manageSliders);

let ang = document.getElementById("angle");
ang.addEventListener("change", manageSliders);
ang.addEventListener("input", manageSliders);
