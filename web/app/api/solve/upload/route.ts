import { NextRequest, NextResponse } from "next/server";
import getConfig from "next/config";

const { serverRuntimeConfig } = getConfig();
const API = serverRuntimeConfig.railwayApiUrl;

export async function POST(req: NextRequest) {
  const method = req.nextUrl.searchParams.get("method") ?? "auto";
  try {
    const body = await req.arrayBuffer();
    const ct = req.headers.get("content-type") ?? "multipart/form-data";
    const upstream = await fetch(
      `${API}/solve/upload?method=${method}`,
      {
        method: "POST",
        headers: { "content-type": ct },
        body,
        signal: AbortSignal.timeout(130_000),
      }
    );
    const data = await upstream.json();
    return NextResponse.json(data, { status: upstream.status });
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : "Upstream error";
    return NextResponse.json({ detail: message }, { status: 502 });
  }
}
