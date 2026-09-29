import { NextRequest, NextResponse } from "next/server";
import getConfig from "next/config";

const { serverRuntimeConfig } = getConfig();
const API = serverRuntimeConfig.railwayApiUrl;

export async function GET(
  _req: NextRequest,
  { params }: { params: { id: string } }
) {
  return NextResponse.json({ error: "Use POST" }, { status: 405 });
}

export async function POST(
  req: NextRequest,
  { params }: { params: { id: string } }
) {
  const method = req.nextUrl.searchParams.get("method") ?? "auto";
  try {
    const upstream = await fetch(
      `${API}/solve/preset/${params.id}?method=${method}`,
      { method: "POST", signal: AbortSignal.timeout(130_000) }
    );
    const data = await upstream.json();
    return NextResponse.json(data, { status: upstream.status });
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : "Upstream error";
    return NextResponse.json({ detail: message }, { status: 502 });
  }
}
